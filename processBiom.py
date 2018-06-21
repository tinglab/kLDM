# -*- coding:utf-8 -*-
from biom import load_table
import os
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description="Filter OTU table and meta data:")
parser.add_argument("-b", "--biom", type=str, required=True, help="The biom file (e.g. ./test/test.biom)")
parser.add_argument("-o", "--result_dir", type=str, required=True,
                    help="Specify the output dir path of the result file")
parser.add_argument("-m", "--meta", type=str, required=True, help="The meta table (samples (rows) * meta (columns)),\
                          corresponding to the sample list file; the first column is sample names and the first row is:\
                          sample_ids    meta1   meta2   ...")

otu_group = parser.add_argument_group("OTUs", description="Filter OTUs")
otu_group.add_argument("-omin", "--otu_min", type=float, default=0.8, \
                       help="Step 1: In every sample, OTUs of which sizes are more than the quantile will be saved (OTUs\
                        with size <= 1 are removed at first)")
otu_group.add_argument("-otop", "--otu_top", type=float, default=0.8, \
                       help="Step 2: OTUs of which average sizes in non-zero samples are more than the quantile\
                        will be saved")
otu_group.add_argument("-not0", "--otu_not_zero", type=float, default=0.25, \
                       help="Step 3: OTUs of which non-zero times in all samples more than not0 * N (all samples number)\
                        will be saved")

sample_group = parser.add_argument_group("Samples", description="Filter samples")
sample_group.add_argument("-smin", "--sample_min", type=float, default=0.10, \
                          help="Step 4.1: Samples of which sizes are less than the quantile will be removed")
sample_group.add_argument("-smax", "--sample_max", type=float, default=0.01, \
                          help="Step 4.2: Samples of which sizes are more than the quantile will be removed")
sample_group.add_argument("-even", "--evenness", type=float, default = 2.0, \
                          help="Step 4.3: Samples of which evenness are less than the value will be removed")

args = parser.parse_args()

print "The biom file is %s" % (args.biom)
print "Filtration Begin~"


class ProcessBiom:
    "Filter samples and OTUs within Biom file"

    def __init__(self, biom_file, meta_file, otu_min, otu_top, otu_not_zero, sample_min, sample_max, result_dir, evenness):
        # load biom file
        try:
            self.table = load_table(biom_file)
        except IOError:
            print "Load biom file failure! Please specify correct biom file and try again~"
            sys.exit()

        print "Load biom file success~"

        self.meta_file = meta_file
        self.otu_top = otu_top
        self.otu_min = otu_min
        self.otu_not_zero = otu_not_zero
        self.sample_min = sample_min
        self.sample_max = sample_max
        self.result_dir = result_dir
        self.evenness = evenness

        # read meta file
        try:
            meta = open(meta_file, 'r')
            col_label = [label.strip().strip("\"") for label in meta.readline().split("\t")]
            col_num = len(col_label)

            sample_name = []
            meta_txt = []
            line = meta.readline()
            while line and len(line) > 0:
                per_line = [per.rstrip() for per in line.split("\t")]
                if len(per_line) != col_num:
                    print "Content Error in meta file! the number of items in some line is not %d. Please check the meta\
                            file and try again~" % (col_num)
                    sys.exit()

                sample_name.append(per_line[0].strip().strip("\""))
                meta_txt.append([float(x) for x in per_line[1:]])
                line = meta.readline()

            self.meta_data = np.array(meta_txt)
            self.meta_names = col_label[1:]
            self.samples_need = sample_name

            print "%d samples are found in meta file" % (len(sample_name))
            print "%d metas are found in meta file: " % (col_num - 1)
            print self.meta_names

        except IOError:
            print "Read meta file failure! Please specify correct meta file and try again~"
            sys.exit()

    # select samples mathced with the sample names in meta data
    def select_samples(self):
        # index of meta data which don't exist in otu table
        remove_index = []
        # samples exist in otu table
        save_index = []
        save_samples = []
        otu_samples = self.table.ids(axis="sample")

        for per in range(len(self.samples_need)):
            per_sample = self.samples_need[per]
            if per_sample in otu_samples:
                save_samples.append(per_sample)
            else:
                remove_index.append(per)

        print "%d samples' ids exist in original OTU table" % (len(save_samples))

        self.table.filter(save_samples, axis="sample")
        samples_select = self.table.ids(axis="sample")

        for per in samples_select:
            index = self.samples_need.index(per)
            save_index.append(index)

        self.meta_data = self.meta_data[np.array(save_index), :]
        self.samples_need = samples_select.tolist()

    # filter otus of which sizes are small and non-zero times are small
    def filter_otus(self):
        (p, n) = self.table.shape
        otu_ids_set = []
        otu_names = np.array(self.table.ids(axis="observation"))
        # Step 1: In every sample, OTUs of which sizes are more than the quantile will be saved
        # 第一步:每个样本中根据OTU大小挑选OTU,过滤掉大小不超过1,且在分位数otu_min一下的OTU
        for per in range(n):
            per_data = self.table[:, per]
            (per_nzero_row, per_nzero_col) = per_data.nonzero()
            per_nzero_data = np.array(per_data[(per_nzero_row, per_nzero_col)].tolist()[0])
            per_mone_index = which(per_nzero_data, lambda x: x > 1)
            if len(per_mone_index) > 0:
                threshold = np.percentile(per_nzero_data[per_mone_index], self.otu_min * 100)
                otu_per_set = which(per_nzero_data, lambda x: x > threshold)
                otu_per_set = per_nzero_row[otu_per_set]
                otu_ids_set.extend(otu_per_set)

        otu_ids_set = set(otu_ids_set)
        otu_ids = otu_names[list(otu_ids_set)]
        self.table.filter(otu_ids, axis="observation")
        otu_ids = self.table.ids(axis="observation")
        print "Step 1: %d OTUs are filtered! %d OTUs are saved!" % (p - len(otu_ids), len(otu_ids))

        # Step 2: OTUs of which average sizes in non-zero samples are more than the quantile will be saved
        # 第二步: 根据OTU在非零样本中的平均大小来过滤,保留大小超过分位数otu_top的OTU
        otu_size = self.table.sum(axis="observation")
        otu_not0 = self.table.nonzero_counts(axis="observation", binary=True)
        otu_size = otu_size / otu_not0
        threshold = np.percentile(otu_size, self.otu_top * 100)
        otu_ids2_index = which(otu_size, lambda x: x > threshold)
        otu_ids2 = otu_ids[otu_ids2_index]
        self.table.filter(otu_ids2, axis="observation")
        otu_ids2 = self.table.ids(axis="observation")
        print "Step 2: %d OTUs are filtered! %d OTUs are saved!" % (len(otu_ids) - len(otu_ids2), len(otu_ids2))

        # Step 3: OTUs of which non-zero times in all samples more than the otu_not_zero * N (all samples number) will be saved
        # 第三步: 根据OTU的非零次数来过滤,过滤掉非零样本数不超过总样本比例的otu_not_zero 的 OTU
        otu_not_zero = self.table.nonzero_counts(axis="observation", binary=True)
        threshold = self.otu_not_zero * n

        otu_ids3_index = which(otu_not_zero, lambda x: x > threshold)
        otu_ids3 = otu_ids2[otu_ids3_index]
        self.table.filter(otu_ids3, axis="observation")
        print "Step 3: %d OTUs are filtered!" % (len(otu_ids2) - len(otu_ids3))
        print "Finally %d OTUs are saved~" % (len(otu_ids3))

    # remove samples with too small sizes or too big sizes
    def filter_samples(self):
        sample_size = self.table.sum(axis="sample")
        # Step 4.1: Samples of which sizes are less than the quantile will be removed
        # 第4.1步 根据样本大小过滤掉size不超过分位数sample_min的样本
        threshold1 = np.percentile(sample_size, self.sample_min * 100)
        # Step 4.2: Samples of which sizes are more than the quantile will be removed
        # 第4.2步 根据样本大小过滤掉size超过分位数sample_max的样本
        threshold2 = np.percentile(sample_size, (1 - self.sample_max) * 100)
        sample_index = which(sample_size, lambda x: threshold1 <= x <= threshold2)
        # Step 4.3: Samples of which evenness are less than the threshold will be removed
        # 第4.3步 根据样本中微生物的相对丰度的均衡度过滤样本
        sample_index_remain = []
        for per_s in sample_index:
            per_sample = self.table.data(self.samples_need[per_s], axis='sample', dense=True)
            per_sample_r = (per_sample + 1) / sum(per_sample + 1)
            if cal_evenness(per_sample_r) > self.evenness:
                sample_index_remain.append(per_s)

        sample_final = np.array(self.samples_need)[sample_index_remain]

        self.table.filter(sample_final, axis="sample")
        samples_final = self.table.ids(axis="sample")
        sample_index = []
        for per in samples_final:
            index = self.samples_need.index(per)
            sample_index.append(index)

        self.meta_data = self.meta_data[np.array(sample_index), :]
        self.samples_final = samples_final

        print "%d samples are saved~" % (len(sample_final))

    # save filtered results into result_dir
    # X N * P
    # M N * Q
    # meta_names Q * 1
    # otu_annotations P * 1
    # sample_names N * 1
    def save_result(self):
        # otu table N * P
        X = self.table.matrix_data.todense().transpose()
        otu_size = np.sum(X, axis=0)
        sort_pair = zip(otu_size.tolist()[0], range(len(otu_size.tolist()[0])))
        sort_pair = sorted(sort_pair, key=lambda x: x[0], reverse=True)
        sort_index = [per[1] for per in sort_pair]

        X = X[:, sort_index]

        # meta table N * Q
        M = self.meta_data
        meta_names = self.meta_names
        otu_ids = self.table.ids(axis="observation")
        otu_names = [";".join(self.table.metadata(axis="observation", id=per)["taxonomy"]) for per in otu_ids]
        otu_names = np.array(otu_names)[sort_index].tolist()

        try:
            if not os.path.exists(self.result_dir):
                os.makedirs(self.result_dir)
            # save the otu table
            np.savetxt(os.path.join(self.result_dir, "otu_table"), X, fmt="%d")
            # save the meta table
            np.savetxt(os.path.join(self.result_dir, "meta_table"), M, fmt="%.8f")
            # save the sample names
            sample_names = open(os.path.join(self.result_dir, "sample_names.txt"), 'w')
            sample_names.write("\n".join(self.samples_final))
            sample_names.close()
            # save the meta names of meta data
            meta_annotation = open(os.path.join(self.result_dir, "meta_annotation.txt"), 'w')
            meta_annotation.write("\n".join(meta_names))
            meta_annotation.close()
            # save the annotations of OTUs
            otu_annotation = open(os.path.join(self.result_dir, "otu_annotation.txt"), 'w')
            otu_annotation.write("\n".join(otu_names))
            otu_annotation.close()
        except IOError:
            print "Save results error!"
            sys.exit()


# return the index in x which satisfy the condition
# the condition is function with bool as return value
def which(x, condition):
    flag = map(condition, x)
    return filter(lambda per: flag[per], range(len(flag)))

def cal_evenness(x):
    return - sum(np.log(x) * x)

processor = ProcessBiom(args.biom, args.meta, args.otu_min, args.otu_top, args.otu_not_zero, args.sample_min, \
                        args.sample_max, args.result_dir, args.evenness)

processor.select_samples()
processor.filter_otus()
processor.filter_samples()
processor.save_result()

print "Filtration Success~"
