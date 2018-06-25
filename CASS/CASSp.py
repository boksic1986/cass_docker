# -*-coding:utf-8-*-
# ! /user/bin/env python
# Author: jiucheng
# Date: 2018/6/5

import json
import logging
import os
import sys
import argparse
from functools import wraps

import sh


logger_name = 'cass'

def logger(level=logging.INFO, name=None, message=None, log_file=None):
    def decorate(func):
        log_name = name if name else func.__module__
        log = logging.getLogger(log_name)
        log.setLevel(logging.INFO)
        if not log.handlers:
            fn = logging.StreamHandler(sys.stderr)
            formatter = logging.Formatter('%(asctime)s %(name)-6s: %(levelname)-4s %(message)s', "%Y-%m-%d %H:%M:%S")
            fn.setFormatter(formatter)
            fn.setLevel(logging.INFO)
            log.addHandler(fn)
        log_msg = message if message else func.__name__

        @wraps(func)
        def wrapper(*args, **kwargs):
            log.log(level=level, msg=log_msg)
            try:
                f = func(*args, **kwargs)
            except sh.ErrorReturnCode :
                raise ValueError ('{}, {}, errors.'.format(log_name, func.__name__))
            else:
                log.log(level=level, msg='{} done.'.format(func.__name__))
                return f
        return wrapper
    return decorate


def lazy_property(func):
    """
    lazy evaluation
    :param func: function
    """
    name = '__lazy__' + func.__name__

    @property
    def lazy(self):
        if hasattr(self, name):
            return getattr(self, name)
        else:
            value = func(self)
            setattr(self, name, value)
            return value

    return lazy


class CassPipe(object):
    """
    cass pipeline in python class
    """

    def __init__(self, input_path, output_path, config_path=None):
        self._input = input_path
        self._out = output_path
        self._config = config_path
        if config_path:
            os.environ['PATH'] += os.pathsep + self.config_path['bio_env']
        self.name = os.path.basename(self._input).split('.')[0]
        self._pipe = self.config_path['pipe_path'] if config_path else os.getenv('PIPE_PATH')
        self.perl = sh.Command(self.config_path['perl']) if config_path else sh.Command('/usr/bin/perl')
        self.R = self.config_path['R'] if config_path else '/usr/bin/R'
        self.convert = sh.Command(self.config_path['convert']) if config_path else sh.Command('/usr/bin/convert')
        self.st = sh.Command('samtools')
        self.bwa = sh.Command('bwa')

    @property
    def pix(self):
        if os.path.exists(self._out):
            return os.path.join(self._out, self.name)
        else:
            os.mkdir(self._out)
            return os.path.join(self._out, self.name)

    @property
    def config_path(self):
        with open(self._config) as cf:
            js = json.load(cf)
        return js


    def pipepl(self, pl):
        return os.path.join(self._pipe, pl)

    def bwa_aln(self):
        bwa = sh.Command('bwa')
        pass

    @logger(name=logger_name, message='start to check bam index')
    def ck_input(self):
        if not os.path.exists(self._input + '.bai'):
            cmd = ['index', '-b', self._input]
            self.st(cmd)

    @logger(name=logger_name, message='start to remove duplicates')
    def drop_dup(self):
        cmd = ['rmdup', self._input, self.pix + '.rmdup.bam']
        self.st(cmd, _out=sys.stdout, _err=sys.stdout)

    @logger(name=logger_name, message='start to sort rmdupped bam files')
    def ca_cnv(self):
        cmd = [self.pipepl('ext.pl'),
               '-i', self.pix + '.rmdup.bam',
               '-o', self.pix + '.ext']
        self.perl(cmd, _err=sys.stdout)

    @lazy_property
    def gender(self):
        cmd = [self.pipepl('tags_gc_gender_info.pl'), '-w',
               self.pipepl('windows_hg19_bwa.aln_500k_slide100k.txt.gz'),
               '-i', self.pix + '.ext.gz',
               '-o', self.pix + '.tagsinfo.gz']
        p = self.perl(cmd, _err=sys.stderr)
        return p.stdout

    @logger(name=logger_name, message='cnv calling for PGS')
    def segment(self):
        cmd = [self.pipepl('segmentation.pl'),
               '-n', self.pipepl('N_region_hg19'),
               '-band', self.pipepl('hg19_cytoBand.txt'),
               '-i', self.pix + '.tagsinfo.gz',
               '-o', self.pix + '.copyratio.cnv.txt',
               '-g', self.gender]
        self.perl(cmd)

    @logger(name=logger_name, message='drop cnv picture')
    def draw_pic(self):
        cmd = [self.pipepl('CNV_drawer.pl'),
               '-r', self.R,
               '-i', self.pix + '.copyratio.cnv.txt',
               '-o', self.pix + '.copyratio',
               '-z', self.pix + '.tagsinfo.gz',
               '-ylim', "0, 3",
               '-name', self.name]
        self.perl(cmd)
        cmd1 = [self.pipepl('boxplot.pl'),
                '-name', self.name,
                '-i', self.pix + '.tagsinfo.gz',
                '-o', self.pix + '.chrratio.xls',
                '-g', self.pix + '.chrratio.pdf',
                '-R', self.R]
        self.perl(cmd1)

    @logger(name=logger_name, message='draw chroms svg')
    def svg(self):
        cmd = [self.pipepl('CNVDrawer.pl'),
               '-lt', 5000000,
               '-i', self.pix + '.copyratio.cnv.txt',
               '-o', self.pix + '.cnv.svg',
               '-gender', self.gender]
        self.perl(cmd)

    @logger(name=logger_name, message='convert svg to png')
    def conv(self):
        cmd = [self.pix + '.cnv.svg', self.pix + '.cnv.png']
        self.convert(cmd, _err=sys.stdout)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", metavar="INPUT BAM/FASTQ", dest='input', required=True),
    parser.add_argument("-o", "--output", metavar="OUTPUT PATH", dest='output', required=True),
    parser.add_argument("-j", "--json", metavar="CONFIG FILE", dest='json')
    return parser.parse_args()



def main():
    args = get_parser()
    file_path = args.input
    out_path = args.output
    json_path = args.json
    cass = CassPipe(file_path, out_path, json_path)
    if file_path.endswith('fq') or file_path.endswith('fastq'):
        cass.bwa_aln()
    elif not file_path.endswith('bam'):
        logging.error('wrong input format')
        sys.exit(0)
    cass.ck_input()
    cass.drop_dup()
    cass.ca_cnv()
    cass.segment()
    cass.draw_pic()
    cass.svg()
    cass.conv()



if __name__ == '__main__':
    main()


