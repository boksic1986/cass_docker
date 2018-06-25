# -*-coding:utf-8-*-
# ! /user/bin/env python
# Author: jiucheng
# Date: 2018/6/21

import os
import sys
import sh
import json
import pickle
import subprocess

from bssh_native_app.BaseSpaceNativeApp import BaseSpaceNativeApp


def get_unique_directory(candidate_directory_path):
    count = 0
    directory_path = candidate_directory_path
    while os.path.isdir(directory_path):
        count += 1
        directory_path = os.path.dirname(candidate_directory_path) + '_' + str(count) + '/'
    return directory_path


def cass(bam_file, output_file):
    candidate_output_directory = os.path.join(output_file, os.path.basename(bam_file).split('.')[0] + '/')
    output_directory_path = get_unique_directory(candidate_output_directory)
    os.makedirs(output_directory_path)
    cmd = ['python', os.path.join(os.getenv('PIPE_PATH'), 'CASSp.py'), '-i', bam_file, '-o', output_directory_path]
    exit_code = subprocess.call(cmd)
    return exit_code, output_directory_path


class CASS(BaseSpaceNativeApp):
    def __init__(self, input_directory_path, output_directory_path, scratch_directory_path, log_directory_path,
                 appsession_id=None):
        BaseSpaceNativeApp.__init__(self, input_directory_path, output_directory_path, scratch_directory_path,
                                    log_directory_path, appsession_id)
        self.__input_directory_path = input_directory_path
        self.__log_directory_path = log_directory_path

    def do_work(self, workspace_directory, appsettings, output_builder):
        file_path = os.path.join(self.__log_directory_path, "AppSession.pk")
        file_path1 = os.path.join(self.__log_directory_path, "Builder.pk")
        with open(file_path, 'wb') as f:
            pickle.dump(appsettings, f)
        with open(file_path1, 'wb') as f1:
            pickle.dump(output_builder, f1)


if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    scratch_dir = sys.argv[3]
    logs_dir = sys.argv[4]
    app = CASS(input_dir, output_dir, scratch_dir, logs_dir)
    app.start()
