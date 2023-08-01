#!/usr/bin/env python3
#  0.1j

import os
import os.path as path
import re
import shutil
import matplotlib.pyplot as plt
import collections
import seaborn as sns

from pyERP import ERP
from pyERP import common

import matplotlib
matplotlib.use('Agg')

OUTPUT_CNV_PLOT = False
OUTPUT_MIXSAMPLE_PLOT = False

# RAW_BASE = {'MES-ab': 8, 'Panel-ab': 8, 'MES-bb': 8, 'Panel-bb': 8, 'WES': 10, 'MITO': 0.7, 'CNVseqPLUS': 20, 'WGS': 90, 'PCRNGS': 0.5, 'HLA': 0.5, }
RAW_BASE = {'MES': 8, 'Panel': 8, 'WES': 10, 'MITO': 0.7, 'CNVseqPLUS': 20, 'WGS': 90, 'PCRNGS': 0.5, 'HLA': 0.5, }
# Q30 = {'MES-ab': 0.8, 'Panel-ab': 0.8, 'MES-bb': 0.8, 'Panel-bb': 0.8, 'WES': 0.8, 'MITO': 0.8, 'CNVseqPLUS': 0.8, 'WGS': 0.8, 'PCRNGS': 0.8, 'HLA': 0.8, }
Q30 = {'MES': 0.8, 'Panel': 0.8, 'WES': 0.8, 'MITO': 0.8, 'CNVseqPLUS': 0.8, 'WGS': 0.8, 'PCRNGS': 0.8, 'HLA': 0.8, }

READS = {'CNVseqPLUS': 120, 'WGS': 600}

# DEPTHAVG = {'MES-ab': 200, 'Panel-ab': 200, 'MES-bb': 200, 'Panel-bb': 200, 'WES': 100, 'MITO': 5000, 'PCRNGS': 2000, 'HLA': 200}
DEPTHAVG = {'MES': 200, 'Panel': 200, 'WES': 100, 'MITO': 5000, 'PCRNGS': 2000, 'HLA': 200}

# COVERAGE = {'MES-ab': 97, 'Panel-ab': 97, 'MES-bb': 97, 'Panel-bb': 97, 'WES': 95, 'WGS': 97, }
COVERAGE = {'MES': 97, 'Panel': 97, 'WES': 95, 'WGS': 97, }

# TOTAl_CNV = {'MES-ab': 1000, 'Panel-ab': 1000}
TOTAl_CNV = {'MES': 2500, 'Panel': 2500}

CNVSEQPLUS_100K = {'CNVseqPLUS': 20, 'WGS': 20, }

P_M_CUTOFF = 50
P_M_RATIO = [0.667, 1.5]

DENOVO = 20


def __judge_qc(rid_dict: dict, detail_dict: dict, qc_dict: dict) -> str:


    data_type = rid_dict['DataType']
    born_flag = rid_dict['req_born']

    # 数据量判定
    qc_raw_base_str = ''
    try:
        if data_type == 'MITO':  # 线粒体项目不做数据量判定
            raise Exception

        if qc_dict['raw_base'] < RAW_BASE[data_type]:
            qc_raw_base_str = '数据量不足{}（cutoff：{}） '.format(qc_dict['raw_base'], RAW_BASE[data_type])
    except:
        qc_raw_base_str = ''

    # Q30判定
    q30_str = ''
    try:
        if qc_dict['q30'] < Q30[data_type]:
            q30_str = 'Q30不足 {}（cutoff：{}）'.format(qc_dict['q30'], Q30[data_type])
    except:
        q30_str = ''

    # 平均深度判定
    depth_str = ''
    try:
        if re.search(r'HLA', rid_dict['GrpItem']) is not None:  # HLA分型项目不做深度判定
            raise Exception

        if qc_dict['depthAvg'] < DEPTHAVG[data_type]:
            depth_str = '平均深度不足 {}（cutoff：{}） '.format(qc_dict['depthAvg'], DEPTHAVG[data_type])
    except:
        depth_str = ''

    # 覆盖度判定
    coverage_str = ''
    try:
        if qc_dict['coverage'] < COVERAGE[data_type]:
            coverage_str = '覆盖度不足（cutoff：{}） '.format(COVERAGE[data_type])
    except:
        coverage_str = ''

    # total_CNV 判定
    total_CNV_str = ''
    try:
        if qc_dict['CNV_total'] > TOTAl_CNV[data_type] and born_flag == '产后':
            total_CNV_str = 'CNV过高 {}（cutoff：{}） '.format(qc_dict['CNV_total'], TOTAl_CNV[data_type])
    except:
        total_CNV_str = ''

    # CNVSEQPLUS 判定
    CNVPLUS_str = ''
    try:
        if qc_dict['CNVSeqPlus_100K'] > CNVSEQPLUS_100K[data_type]:
            CNVPLUS_str = '100K片段CNV过高 {}（cutoff：{}） '.format(qc_dict['CNVSeqPlus_100K'], CNVSEQPLUS_100K[data_type])
    except:
        CNVPLUS_str = ''

    # 性别 判定
    sex_str = ''
    try:
        if data_type == 'PCRNGS':  # PCRNGS项目 不做性别判定
            raise Exception

        if qc_dict['sex_record'] != '未知' and qc_dict['sex_calculate'] is not None and qc_dict['sex_record'] != qc_dict['sex_calculate']:
            sex_str = '性别不符 '
    except:
        sex_str = ''

    # 家系 判定
    pedigree_str = ''
    p = qc_dict['P_count']
    m = qc_dict['M_count']
    try:
        if data_type == 'PCRNGS':  # PCRNGS项目不做家系判定
            raise Exception

        if p is None or m is None:
            pedigree_str = ''
        elif p < P_M_CUTOFF or m < P_M_CUTOFF:
            pedigree_str = '父系突变{p}，母系突变{m} (cutoff：{cutoff}) '.format(p = p,
                                                                                          m = m,
                                                                                          cutoff = P_M_CUTOFF)
        elif p / m <= P_M_RATIO[0] or p / m >= P_M_RATIO[1]:
            pedigree_str = '父母比例差距过大 {ratio}（cutoff：{cutoff}） '.format(ratio = round(p / m, 2),
                                                                                 cutoff = P_M_RATIO)
        else:
            raise ValueError
    except:
        pedigree_str = ''


    # deNovo 判定
    deNovo_str = ''
    try:
        if qc_dict['deNovo'] is None:
            deNovo_str = ''
        elif qc_dict['deNovo'] > DENOVO:
            deNovo_str = '（新发突变过高{deNovo} cutoff：{cutoff}） '.format(deNovo = qc_dict['deNovo'],
                                                                               cutoff = DENOVO)
    except:
        deNovo_str = ''



    return sex_str + pedigree_str + qc_raw_base_str + q30_str + depth_str + coverage_str + total_CNV_str + CNVPLUS_str + deNovo_str


def output_qc(file_handle, rid_dict: dict, detail_dict: dict, qc_dict: dict, no: int, print_header = False, custom = '') -> int:

    header = '序号\t质控\t批次号\t流水号\t项目号\t家庭号\t姓名\t家庭关系\t项目\t数据类型\t产前后\t样本类型\t速度\t数据量(G)\tQ30\t捕获率\t平均深度\t覆盖度\t正反向距离\t错误率\tTotal_CNV\t整倍体数\t三体染色体\tCNVSeq_totalSample\tCNVSeq_totalAbnormal\tCNVSeqPlus_100K\tCNVSeqPlus_Trip\t录入性别\t计算性别\t核型性别\tP_count\tM_count\tPvsM\tdeNovo\t病历'
    if print_header:
        file_handle.writelines(header + '\n')

    if rid_dict is None and detail_dict is None and qc_dict is None:  # 如果无法获取质控信息
        line_str = '{no}\t{custom}'.format(no = no, custom = custom.strip())
    else:
        try:
            qc_str = __judge_qc(rid_dict, detail_dict, qc_dict)
            line_str = '{no}\t{qc_str}\t{BatchNo}\t{BarCode}\t{GrpEncode}\t{Family_id}\t{PatName}\t{Family_relation}\t{GrpItem}\t{DataType}\t{Born_flag}\t{sp_name}\t{urgent}\t{Raw_bases}\t{Q30:1.1%}\t{On_target_rate:1.1%}\t{DepthAvg}±{DepthSD}\t{Coverage}%\t{Insert_size_average}\t{error_rate:1.2%}\t{CNV_total}\t{Ploidy}\t{Trisome}\t{CNVSeq_totalSample}\t{CNVSeq_totalAbnormal}\t{CNVSeqPlus_100K}\t{CNVSeqPlus_Trip}\t{Sex_record}\t{Sex_calculate}\t{Sex_info}\t{P_count}\t{M_count}\t{PvsM}\t{deNovo}\t{p_clinical_sign}\t{custom}'.format(
                no = no,
                qc_str = qc_str,
                BatchNo = '-' if rid_dict['BatchNo'] is None else rid_dict['BatchNo'],
                BarCode = rid_dict['BarCode'],
                GrpEncode = rid_dict['GrpEncode'],
                Family_id = rid_dict['p_family_no'],
                PatName = rid_dict['PatName'].strip(),
                Family_relation = rid_dict['family_relation'],
                GrpItem = rid_dict['GrpItem'].strip(),
                DataType = rid_dict['DataType'],
                Born_flag = rid_dict['req_born'],
                sp_name = rid_dict['sp_name'],
                urgent = '-' if detail_dict['urgent'] is None else detail_dict['urgent'],

                Raw_bases = qc_dict['raw_base'],
                Q30 = -1 if qc_dict['q30'] is None else qc_dict['q30'],
                On_target_rate = -1 if qc_dict['on_target_rate'] is None else qc_dict['on_target_rate'],
                DepthAvg = qc_dict['depthAvg'],
                DepthSD = qc_dict['depthSD'],
                Coverage = qc_dict['coverage'],
                Insert_size_average = qc_dict['insert_size_average'],
                error_rate = -1 if qc_dict['error_rate'] is None else qc_dict['error_rate'],

                CNV_total = '-' if qc_dict['CNV_total'] is None else qc_dict['CNV_total'],
                Ploidy = '-' if qc_dict['ploidy'] is None or qc_dict['ploidy'] == 'Diploid' else qc_dict['ploidy'],
                Trisome = '-' if qc_dict['trisome'] is None or qc_dict['trisome'] == 'None' else qc_dict['trisome'],
                CNVSeq_totalSample = '-' if qc_dict['CNVSeq_totalSample'] is None else qc_dict['CNVSeq_totalSample'],
                CNVSeq_totalAbnormal = '-' if qc_dict['CNVSeq_totalAbnormal'] is None else qc_dict['CNVSeq_totalAbnormal'],
                CNVSeqPlus_100K = '-' if qc_dict['CNVSeqPlus_100K'] is None else qc_dict['CNVSeqPlus_100K'],
                CNVSeqPlus_Trip = '-' if qc_dict['CNVSeqPlus_Trip'] is None else qc_dict['CNVSeqPlus_Trip'],
                Sex_calculate = '-' if qc_dict['sex_calculate'] is None else qc_dict['sex_calculate'],
                Sex_record = '-' if detail_dict['sex_recorded'] is None else detail_dict['sex_recorded'],
                Sex_info = '-' if qc_dict['sexInfo'] is None else qc_dict['sexInfo'],
                P_count = '-' if qc_dict['P_count'] is None else qc_dict['P_count'],
                M_count = '-' if qc_dict['M_count'] is None else qc_dict['M_count'],
                PvsM = '-' if qc_dict['PvsM'] is None else qc_dict['PvsM'],
                deNovo = '-' if qc_dict['deNovo'] is None else qc_dict['deNovo'],
                p_clinical_sign = '-' if detail_dict['p_clinical_sign'] is None else detail_dict['p_clinical_sign'].strip(),
                custom = custom.strip()
            )
        except Exception:
            line_str = str(no)

    # print(' ' * 30, end = '\r')
    # print(no, rid_dict['BarCode'], rid_dict['GrpEncode'])
    file_handle.writelines(line_str + '\n')
    return 0

def output_mixsample_plot(mixsample_lst: list, plot: str) -> int:
    try:
        plt.figure()
        histplot = sns.histplot(data = mixsample_lst, stat = 'probability', binwidth = 1, kde = True).get_figure()
        histplot.savefig(plot)
        plt.close()
    except:
        print('无法生成{}'.format(plot))


    return 0

def main(arguments_dict = {}):

    e = ERP.ERP(username = 'tianran.yao', password = 'qbTgqaH2C6yfJQzz')
    # r = e._api_GetRequestNote(r_code = '1000151877')

    file_in_str = ''
    while file_in_str == '':
        file_in_str = input('输入列表文件名：')
    file_out_str = path.splitext(file_in_str)[0] + '_qc.txt'
    folder_out_str = path.splitext(file_in_str)[0] + '_qc'

    plot = input('是否下载图片（输入c，下载CNV图片。输入m，下载混样图。输入a，下载所有。回车不下载）:')
    if plot == 'c':
        OUTPUT_CNV_PLOT = True
        OUTPUT_MIXSAMPLE_PLOT = False
    elif plot == 'm':
        OUTPUT_CNV_PLOT = False
        OUTPUT_MIXSAMPLE_PLOT = True
    elif plot == 'a':
        OUTPUT_CNV_PLOT = True
        OUTPUT_MIXSAMPLE_PLOT = True
    else:
        OUTPUT_CNV_PLOT = False
        OUTPUT_MIXSAMPLE_PLOT = False


    if os.access(file_out_str, mode = os.R_OK):
        os.remove(file_out_str)

    if not os.access(folder_out_str, mode = os.R_OK):
        os.mkdir(folder_out_str)


    # ==================get qc info =====================
    in_f = open(file_in_str)
    out_f = open(file_out_str, 'at', buffering = 1)
    i = 0
    for line_str in in_f.readlines():
        i += 1

        if line_str == '' or line_str.startswith('#'):
            continue

        try:
            row_lst = re.split(r'[,\s]+', line_str.strip())
            barcode = row_lst[0]
            grpcode = row_lst[1]
        except:
            message = '第{}行 {} 格式错误. 跳过!'.format(i, line_str)
            print(message)
            continue

        print('{} {}-{} ...'.format(i, barcode, grpcode), end = '\r')
        # 质控查询 rid_dict
        try:
            index = '{}-{}'.format(barcode, grpcode)
            rid_dict = e._api_SearchByPage(bar_code = barcode, grp_encoded = grpcode)
            rid_dict = rid_dict[index]
            rid = rid_dict['rid']
            data_type = rid_dict['DataType']
        except:
            rid_dict = collections.defaultdict(lambda: None)
            message = '{} 无质控信息            '.format(index)
            print(message)
            if i == 1:
                r = output_qc(out_f, None, None, None, i, print_header = True, custom = message)
            else:
                r = output_qc(out_f, None, None, None, i, print_header = False, custom = message)
            continue

        # 质控详情 qc_dict
        e.reset()
        try:
            qc_dict = e._api_GetQcInfo(rid)
        except:
            qc_dict = collections.defaultdict(lambda: None)
            message = '{} {} 无质控详情            '.format(i, index)
            print(message)

        # 申请单详情页面  detail_dict
        try:
            detail_dict = e._api_GetRequestNote(rid = rid)  # extract detail
        except:
            detail_dict = collections.defaultdict(lambda: None)
            message = '{} {} 无申请单详情            '.format(i, index)
            print(message)

        # 获取CNV和UDP图
        if OUTPUT_CNV_PLOT:
            try:
                file_lst = e._get_CNV_UPD_plot(barcode, grpcode, folder_out_str, prefix = 'CNV-' + str(i).rjust(3, '0'))
                if len(file_lst) == 0:
                    raise ConnectionError
            except:
                message = '{} {} 无法获取CNV或UPD图            '.format(i, index)
                print(message)

        # 质控详情（混样信息）
        if OUTPUT_MIXSAMPLE_PLOT:
            try:
                data_id = rid_dict['data_id']
                if data_type in ['MES', 'Panel', 'MITO', 'WES', 'WGS']:
                    mixsample_lst = e._api_GetMixSampleQc(data_id)
                    r = output_mixsample_plot(mixsample_lst, '{folder}/mixsample_{i:03}_{index}.png'.format(i = i,
                                                                                                            folder = folder_out_str,
                                                                                                            index = index)
                                              )
            except:
                message = '{} {} 无混样信息            '.format(i, index)
                print(message)

        if i == 1:
            r = output_qc(out_f, rid_dict, detail_dict, qc_dict, i, print_header = True)
        else:
            r = output_qc(out_f, rid_dict, detail_dict, qc_dict, i, print_header = False)
    in_f.close()
    out_f.close()

    try:
        shutil.copy(file_out_str, path.join(folder_out_str, file_out_str))
    except:
        print(file_out_str, path.join(folder_out_str, file_out_str))

    print('完成')

    return

if __name__ == '__main__':
    # arguments_dict = get_arguments()
    main()

