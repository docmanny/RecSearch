""" Common system for reading and writing BED file format"""

from pathlib import Path
import re

def read(file_location):
    bedfile = Path(file_location)
    assert bedfile.exists(), "Given bedfile path does not exist!"
    assert bedfile.is_file(), "Given bedfile path was not a file! Did you provide a directory?"
    field_names = ["chr", "start", "end", "name",
                   "score", "strand", "thick_start", "thick_end",
                   "item_rgb", "block_count", "block_sizes", "block_starts",
                   "detail_id", "detail_description"]
    header_names = ['name', 'description', 'type', 'visibility', 'color', 'itemRgb', 'colorByStrand', 'useScore',
                    'group', 'priority', 'db', 'offset', 'maxItems', 'url', 'htmlUrl', 'bigDataUrl']
    pats = [r'{0}="?[A-Za-z0-9 ]+"?(?=\s)'.format(i) for i in header_names]
    pats[-3:] = [r'{0}=".+"(?=\s)'.format(i) for i in header_names[-3:]]
    pats_re = [re.compile(p) for p in pats]
    header_dict = {}
    is_bed_detail=False
    with bedfile.open() as bed:
        for line in bed:
            #print(line)
            if line.startswith("browser"):
                continue
            elif line.startswith("track"):
                items = list(filter(None, [p.findall(line) for p in pats_re]))
                p2 = re.compile("(\w+)=\"?(.+)\"?")
                header_dict = {k: v for k, v in [p2.match(i[0]).groups() for i in items]}
                #for k, v in header_dict.items():
                    #print(k,v, sep="\t")
                if "type" in header_dict:
                    is_bed_detail = header_dict["type"] == "bedDetail"
                break
            else:
                break
        records = (line.strip().split("\t") for line in bed)
        if is_bed_detail:
            """
            record = next(records)
            core_record = record[:-2]
            details = record[-2:]
            rec_dict = {k:v for k, v in zip(field_names, core_record)}
            rec_dict.update({k:v for k, v in zip(field_names[-2:], details)})
            for k, r in rec_dict.items():
                print(k, r, sep="\t")
            bl = BEDLine(**rec_dict)
            print(bl)
            bed_record = BEDRecord(bl)
            for record in records:
                core_record = record[:-2]
                details = record[-2:]
                rec_dict = {k: v for k, v in zip(field_names, core_record)}
                rec_dict.update({k: v for k, v in zip(field_names[-2:], details)})
                bed_record = bed_record + BEDLine(**rec_dict)"""
            raise NotImplementedError("The BED Detail format is a mess and needs clarification. Right now, the ad"
                                      "verbatim description of what it entails according to UCSC does not match their"
                                      "own example, so until further notice, there will be no reading in a BED Detail"
                                      "file directly.")
        else:
            bed_record = BEDRecord(*[BEDLine(**{k: v for k, v in zip(field_names, record)}) for record in records])
        if header_dict:
            for key in header_dict:
                bed_record.header[key] = header_dict[key]
        else:
            bed_record.header["name"] = bedfile.name
    return bed_record




class BEDLine(object):
    def __init__(self, **kwargs):
        self.chr = kwargs.pop("chr", ".")
        self.start = int(kwargs.pop("start", 0))
        self.end = int(kwargs.pop("end", 1))
        self.name = kwargs.pop("name", ".")
        self.score = int(kwargs.pop("score", 0))
        self.strand = kwargs.pop("strand", ".")
        self.thick_start = int(kwargs.pop("thick_start", self.start))
        self.thick_end = int(kwargs.pop("thick_end", self.end))
        self.item_rgb = kwargs.pop("item_rgb", "0,0,0")
        self.block_count = int(kwargs.pop("block_count", 0))
        self.block_sizes = kwargs.pop("block_sizes", self.end-self.start)
        self.block_starts = kwargs.pop("block_starts", self.start)
        if "detail_id" in kwargs or "detail_description" in kwargs:
            self.detail_id = kwargs.pop("detail_id", ".")
            self.detail_description = kwargs.pop("detail_description", ".")

    def is_bed_detail(self):
        return all([hasattr(self, attr) for attr in ("detail_id", "detail_description")])

    def is_bed12(self):
        return all([hasattr(self, attr) for attr in ["block_starts", "block_sizes", "block_count", "item_rgb",
                                                     "thick_end", "thick_start", "strand", "score",
                                                     "name", "end", "start","chr"]])

    def is_bed6(self):
        return all([hasattr(self, attr) for attr in ["strand", "score", "name", "end", "start", "chr"]])

    def is_bed4(self):
        return all([hasattr(self, attr) for attr in ["name", "end", "start", "chr"]])

    def as_bed_detail(self, detail_id, detail_description):
        self.detail_id = detail_id
        self.detail_description = detail_description

    def as_bed12(self):
        for attr in ["block_starts", "block_sizes", "block_count", "item_rgb",
                     "thick_end", "thick_start", "strand", "score",
                     "name", "end", "start","chr"]:
            if hasattr(self, attr):
                continue
            else:
                self.__setattr__(attr, ".")

    def as_bed6(self):
        for attr in ["block_starts", "block_sizes", "block_count",
                     "item_rgb", "thick_end", "thick_start"]:
            try:
                self.__delattr__(attr)
            except AttributeError:
                continue
        for attr in ["chr", "start", "end", "name", "score", "strand"]:
            if hasattr(self, attr):
                continue
            else:
                self.__setattr__(attr, ".")


    def as_bed4(self):
        for attr in ["block_starts", "block_sizes", "block_count", "item_rgb",
                     "thick_end", "thick_start", "strand", "score"]:
            try:
                self.__delattr__(attr)
            except AttributeError:
                continue
        for attr in ["chr", "start", "end", "name"]:
            if hasattr(self, attr):
                continue
            else:
                self.__setattr__(attr, ".")

    def bed_type(self):
        if self.is_bed12():
            return "BED12"
        elif self.is_bed6():
            return "BED6"
        elif self.is_bed4():
            return "BED4"
        else:
            return "BED"

    def __str__(self):
        bed = [str(i) for i in [self.chr, self.start, self.end, self.name, self.score, self.strand, self.thick_start,
                                self.thick_end, self.item_rgb, self.block_count, self.block_sizes, self.block_starts]]
        if self.is_bed_detail():
            bed.append(str(self.detail_id))
            bed.append(str(self.detail_description))
        return "\t".join(bed)

    def __repr__(self):
        r = "BEDLine("
        r += ", ".join("{0}={1}".format(k, str(v)) for k,v in self.__dict__.items())
        r += ")"
        return r


class BEDRecord(list):
    def __init__(self, *bed_lines):
        super(list, self).__init__()
        self.header = dict(name=None, description=None, type=None, visibility=None, color=None, itemRgb=None,
                           colorByStrand=None, useScore=None, group=None, priority=None, db=None, offset=None,
                           maxItems=None, url=None, htmlUrl=None, bigDataUrl=None)
        if isinstance(bed_lines, BEDLine):
            self.bed_type = bed_lines.bed_type()
            self.__add__(bed_lines)
        else:
            for bed_line in bed_lines:
                assert isinstance(bed_line, BEDLine), "Only BEDLine objects can be added to a BED Record!"
                self.bed_type = bed_line.bed_type()
                self.__add__(bed_line)

    def is_bed_detail(self):
        return self.header["type"] == "bedDetail"

    def as_bed4(self):
        self.bed_type = "BED4"
        for line in self:
            line.as_bed4()

    def as_bed6(self):
        self.bed_type = "BED6"
        for line in self:
            line.as_bed6()

    def as_bed12(self):
        self.bed_type = "BED12"
        for line in self:
            line.as_bed12()

    def as_bed_detail(self):
        self.header["type"] = "bedDetail"
        for line in self:
            line.as_bed_detail(".", ".")

    def ids(self):
        return (line.name for line in self)

    def write(self, file_location, with_header=False):
        bedfile = Path(file_location)
        with bedfile.open("w") as outf:
            if with_header:
                outf.write(("track"
                            " ".join(('{0}="{1}"'.format(k,v) for k,v in self.header.items()))
                            ))
                outf.write("\n")
            for line in self:
                outf.write(str(line))
                outf.write("\n")

    def sort(self, key=False, reverse=False):
        if key:
            super(BEDRecord, self).sort(key=key, reverse=reverse)
        else:
            super(BEDRecord, self).sort(key=lambda line: line.name, reverse=reverse)
            super(BEDRecord, self).sort(key=lambda line: line.end, reverse=reverse)
            super(BEDRecord, self).sort(key=lambda line: line.start, reverse=reverse)
            super(BEDRecord, self).sort(key=lambda line: line.chr, reverse=reverse)

    def __add__(self, bed_line):
        assert isinstance(bed_line, BEDLine), "Only BEDLine objects can be added to a BED Record!"
        if self.is_bed_detail() != bed_line.is_bed_detail():
            if self.is_bed_detail():
                bed_line.as_bed_detail(".", ".")
            else:
                self.as_bed_detail()
        if self.bed_type == bed_line.bed_type():
            return self.append(bed_line)
        else:
            try:
                if self.bed_type == "BED4":
                    if bed_line.bed_type() == "BED6":
                        self.as_bed6()
                        self.append(bed_line)
                    elif bed_line.bed_type() == "BED12":
                        self.as_bed12()
                        self.append(bed_line)
                    else:
                        raise AssertionError
                elif self.bed_type == "BED6":
                    if bed_line.bed_type() == "BED4":
                        bed_line.as_bed4()
                        self.append(bed_line)
                    elif bed_line.bed_type() == "BED12":
                        self.as_bed12()
                        self.append(bed_line)
                    else:
                        raise AssertionError
                elif self.bed_type == "BED12":
                    if bed_line.bed_type() == "BED4":
                        bed_line.as_bed12()
                        self.append(bed_line)
                    elif bed_line.bed_type() == "BED6":
                        bed_line.as_bed12()
                        self.append(bed_line)
                    else:
                        raise AssertionError
                else:
                    raise AssertionError
            except AssertionError:
                return AssertionError("Cannot add a BEDLine object of type {0} "
                                      "to a BEDRecord object of type {1}!".format(bed_line.bed_type(),
                                                                                  self.bed_type))

    def __radd__(self, other):
        self.__add__(other)

    def __repr__(self):
        return "BedRecord(type={0}, length={1})".format(self.bed_type, len(self))

    def __str__(self):
        return "\n".join(str(i) for i in self)
