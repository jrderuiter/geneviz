from __future__ import absolute_import, division, print_function

#pylint: disable=W0622,W0611
from builtins import (ascii, bytes, chr, dict, filter, hex, input, int, map,
                      next, oct, open, pow, range, round, str, super, zip)
from future.utils import native_str

import contextlib
import itertools
import os
import subprocess

import pysam
import numpy as np
import pandas as pd


class TabixFile(object):
    def __init__(self, filename, parser, **kwargs):
        self._filename = filename
        self._parser = parser
        self._kws = kwargs

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              filters=None,
              incl_left=True,
              incl_right=True,
              raw=False):

        # Fetch and filter records.
        records = self._fetch_records(reference, start, end)
        records = self._filter_records(records, filters, start, end, incl_left,
                                       incl_right)

        # Yield records.
        if not raw:
            records = (self._to_series(rec) for rec in records)

        for record in records:
            yield record

    def _fetch_records(self, reference=None, start=None, end=None, **kwargs):
        # Open file.
        file_obj = pysam.TabixFile(
            native_str(self._filename), parser=self._parser, **self._kws)

        with contextlib.closing(file_obj) as file_:
            # Ensure reference is a native str.
            if reference is not None:
                reference = native_str(reference)

            # For some reason pysam does not fetch all records if reference
            # is None under Python 2.7. To fix this, here we simply chain all
            # the contig records into one iterable.
            if reference is None:
                records = (file_.fetch(
                    reference=contig, **kwargs) for contig in file_.contigs)
                records = itertools.chain.from_iterable(records)
            else:
                records = file_.fetch(
                    reference=reference, start=start, end=end, **kwargs)

            for record in records:
                yield record

    def _filter_records(self,
                        records,
                        filters=None,
                        start=None,
                        end=None,
                        incl_left=True,
                        incl_right=True):

        # Filter records on additional filter values.
        if filters is not None:

            def _filter(records, name, value):
                for r in records:
                    if hasattr(r, name) and getattr(r, name) == value:
                        yield r

            for name, value in filters.items():
                records = _filter(records, name, value)

        # Filter inclusive/exclusive ranges if needed.
        if not incl_left:
            records = filter(lambda r: r.start > start, records)

        if not incl_right:
            records = filter(lambda r: r.end < end, records)

        for record in records:
            yield record

    def fetch_frame(self,
                    reference=None,
                    start=None,
                    end=None,
                    filters=None,
                    incl_left=True,
                    incl_right=True):
        records = self.fetch(
            reference,
            start,
            end,
            filters=filters,
            incl_left=incl_left,
            incl_right=incl_right)
        return self._frame_constructor().from_records(records)

    @classmethod
    def _to_series(cls, record):
        raise NotImplementedError()

    @classmethod
    def _frame_constructor(cls):
        raise NotImplementedError()


class TabixFrame(pd.DataFrame):
    @property
    def _constructor(self):
        raise NotImplementedError()

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              filters=None,
              incl_left=True,
              incl_right=True,
              raw=False):

        subset = self.fetch_frame(reference, start, end, filters, incl_left,
                                  incl_right)

        if raw:
            records = subset.itertuples()
        else:
            records = (row for _, row in subset.iterrows())

        for rec in records:
            yield rec

    def fetch_frame(self,
                    reference=None,
                    start=None,
                    end=None,
                    filters=None,
                    incl_left=True,
                    incl_right=True):

        return _get_region(
            self,
            reference,
            start,
            end,
            filters=filters,
            incl_left=incl_left,
            incl_right=incl_right)


def _parse_float(value):
    try:
        return float(value)
    except ValueError:
        return np.nan


class GtfFile(TabixFile):

    TYPE_MAP = {3: int, 4: int, 5: _parse_float}

    FIELDS = ('contig', 'source', 'feature', 'start', 'end', 'score', 'strand',
              'frame', 'attribute')

    def __init__(self, filename, **kwargs):
        file_path = str(filename)
        if not file_path.endswith('.gz'):
            if os.path.exists(file_path + '.gz'):
                file_path += '.gz'
            else:
                file_path = self.compress(file_path, create_index=True)

        super().__init__(filename, parser=pysam.asGTF(), **kwargs)

    @classmethod
    def _to_series(cls, record):
        rec_values = tuple((cls.TYPE_MAP.get(i, lambda x: x)(val)
                            for i, val in enumerate(record)))
        attr_keys, attr_values = zip(*dict(record).items())
        return pd.Series(
            rec_values[:-1] + attr_values, index=cls.FIELDS[:-1] + attr_keys)

    @classmethod
    def _frame_constructor(cls):
        return GtfFrame

    def get_gene(self,
                 gene_id,
                 feature_type='gene',
                 field_name='gene_id',
                 **kwargs):
        # Add feature filter to filters (if given).
        filters = kwargs.pop('filter', {})
        filters['feature'] = feature_type

        # Search for gene record.
        records = self.fetch(filters=filters, raw=True, **kwargs)
        for record in records:
            if record[native_str(field_name)] == gene_id:
                return self._to_series(record)

        raise ValueError('Gene {} does not exist'.format(gene_id))

    @classmethod
    def compress(cls, filename, out_path=None, sort=True, create_index=True):
        """Compresses and indexes a gtf file using bgzip and tabix."""

        # Base output path on original file name.
        out_path = out_path or filename + '.gz'

        if sort:
            # Sort file before compressing and indexing.
            filename = cls.sort(filename, out_path=out_path + '.tmp')

        # Gzip and index file.
        gzipped_path = bgzip(filename, out_path=out_path)

        if create_index:
            tabix(gzipped_path, preset='gff')

        # Clean up temp file.
        if filename.endswith('.tmp'):
            os.unlink(filename)

        return gzipped_path

    @classmethod
    def sort(cls, filename, out_path):
        """Sorts a gtf file by position, as required for tabix."""
        with open(out_path, 'w') as out_file:
            cmd = '(grep ^"#" {0}; grep -v ^"#" {0} ' '| sort -k1,1 -k4,4n)'
            subprocess.check_call(
                cmd.format(filename), stdout=out_file, shell=True)
        return out_path

    def __repr__(self):
        return '<GtfFile filename={!r}>'.format(self._filename)


class GtfFrame(TabixFrame):
    @property
    def _constructor(self):
        return GtfFrame

    @classmethod
    def read_csv(cls, path, *args, **kwargs):
        frame = pd.read_csv(path, *args, sep='\t', comment='#', **kwargs)
        return cls._format_frame(frame)

    @classmethod
    def from_records(cls, data, *args, **kwargs):
        # Build frame.
        frame = super().from_records(data, *args, **kwargs)

        # Handle empty case.
        if len(frame) == 0:
            frame = pd.DataFrame([], columns=GtfFile.FIELDS[:-1])

        return cls._format_frame(frame)

    @classmethod
    def _format_frame(cls, frame):
        # Convert some columns to categorical.
        frame['contig'] = frame['contig'].astype('category')
        frame['feature'] = frame['feature'].astype('category')
        frame['source'] = frame['source'].astype('category')
        frame['frame'] = frame['frame'].astype('category')

        # Order columns to start with normal attributes.
        frame = _reorder_columns(frame, GtfFile.FIELDS[:-1])

        return frame

    def get_gene(self, gene_id):
        result = self.ix[((self['feature'] == 'gene') &
                          (self['gene_id'] == gene_id))]

        if len(result) == 0:
            raise ValueError('Gene {} does not exist'.format(gene_id))

        return result


class BedFile(TabixFile):
    TYPE_MAP = {1: int, 2: int, 4: _parse_float}

    FIELDS = ('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
              'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
              'blockStarts')

    def __init__(self, filename, **kwargs):
        filename = str(filename)
        if not filename.endswith('.gz'):
            if os.path.exists(filename + '.gz'):
                filename += '.gz'
            else:
                filename = self.compress(filename, create_index=True)

        super().__init__(filename, parser=pysam.asBed(), **kwargs)

    @classmethod
    def _to_series(cls, record):
        values = tuple((cls.TYPE_MAP.get(i, lambda x: x)(val)
                        for i, val in enumerate(record)))
        return pd.Series(values, index=cls.FIELDS[:len(values)])

    @classmethod
    def _frame_constructor(cls):
        return BedFrame

    @classmethod
    def compress(cls, filename, out_path=None, sort=True, create_index=True):
        """Compresses and indexes a bed file using bgzip and tabix."""

        # Base output path on original file name.
        out_path = out_path or filename + '.gz'

        if sort:
            # Read and sort bed file.
            frame = BedFrame.read(filename)
            frame.sort_values(inplace=True)

            filename = out_path + '.tmp'
            frame.write(filename)

        # Gzip and index file.
        gzipped_path = bgzip(filename, out_path=out_path)

        if create_index:
            tabix(gzipped_path, preset='bed')

        # Clean up temp file.
        if filename.endswith('.tmp'):
            os.unlink(filename)

        return gzipped_path

    @classmethod
    def sort(cls, filename, out_path):
        """Sorts a gtf file by position, as required for tabix."""
        with open(out_path, 'w') as out_file:
            cmd = '(grep ^"#" {0}; grep -v ^"#" {0} ' '| sort -k1,1 -k4,4n)'
            subprocess.check_call(
                cmd.format(filename), stdout=out_file, shell=True)
        return out_path


class BedFrame(TabixFrame):

    COL_NAMES = ('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                 'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                 'blockSizes', 'blockStarts')

    @property
    def _constructor(self):
        return BedFrame

    @classmethod
    def read(cls, file_path):
        with open(str(file_path), 'r') as file_:
            skip = 1 if next(iter(file_)).startswith('track') else 0

        return cls(
            pd.read_csv(
                str(file_path),
                sep='\t',
                skiprows=skip,
                names=cls.COL_NAMES,
                dtype={'chrom': str}))

    def write(self, file_path):
        return self.to_csv(
            str(file_path), sep=native_str('\t'), header=False, index=False)

    def sort_values(self, *args, **kwargs):
        if len(args) > 0:
            by = args[0]
        else:
            by = kwargs.get('by', ['chrom', 'chromStart', 'chromEnd'])

        return super().sort_values(by, *args, **kwargs)

    @classmethod
    def from_records(cls, data, *args, **kwargs):
        # Build frame.
        frame = super().from_records(data, *args, **kwargs)

        # Handle empty case.
        if len(frame) == 0:
            frame = pd.DataFrame([], columns=BedFile.FIELDS)

        # Order columns.
        frame = _reorder_columns(frame, BedFile.FIELDS[:frame.shape[1]])

        return frame


def _reorder_columns(frame, order):
    columns = list(order)
    extra_columns = sorted([c for c in frame.columns if c not in set(columns)])
    return frame[columns + extra_columns]


def _get_region(frame,
                reference,
                start=None,
                end=None,
                filters=None,
                incl_left=True,
                incl_right=True,
                ref_col='contig',
                start_col='start',
                end_col='end'):
    # Filter on passed range.
    mask = frame[ref_col] == reference

    if start is not None:
        mask &= frame[start_col] <= end

    if end is not None:
        mask &= frame[end_col] >= start

    # Filter for inclusiveness.
    if not incl_left:
        mask &= frame['start'] > start

    if not incl_right:
        mask &= frame['end'] < end

    # Apply any additional filters.
    if filters is not None:
        for name, value in filters.items():
            mask &= frame[name] == value

    return frame.ix[mask]


def bgzip(file_path, out_path=None):
    if out_path is None:
        out_path = file_path + '.gz'

    with open(out_path, 'w') as out_file:
        subprocess.check_call(['bgzip', '-c', str(file_path)], stdout=out_file)

    return out_path


def tabix(file_path, preset):
    subprocess.check_call(['tabix', '-p', preset, file_path])
