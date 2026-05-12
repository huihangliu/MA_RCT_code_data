import numpy as np
import os
import re

file_path = 'D:/Documents/Codes/R/RCT-R/'

pathDir = os.listdir(file_path)

# Check each file's line number.
for file_name in pathDir:
    if 'copied_from_writting.rmd' not in file_name:
        continue
    print(file_name)

    with open(file_path + file_name, "r", encoding='utf-8') as file_in, open(file_path + 'tmp.Rmd', "w", encoding='utf-8') as file_out:
        file_lines = file_in.readlines()
        for cur_line in file_lines:
        # for num_line in range(0,len(file_lines)):
            # # get current line
            # cur_line = file_lines[num_line]
            text = cur_line

            # change YAML
            if not re.search('layout:[\S\s]*post\n', text) is None:
                text = 'output:\t\tbookdown::html_document2\n'

            # delete one # from titles
            text = text.replace("## ", "# ")
            text = text.replace("### ", "## ")
            text = text.replace("#### ", "### ")
            text = text.replace("##### ", "#### ")

            # in-line equation
            text = text.replace("$$", "$")

            # break-line equation
            text = text.replace("$\n", "\n")  ## I put break-line equation into a align envir. so $$ is not needed. 

            # equation label
            tmp = re.search('\\\\label{[\S\s]*?}', text)
            if not tmp is None:
                text = text[0:tmp.span()[0]] + '(\\#' + text[(tmp.span()[0]+7):(tmp.span()[1]-1)] + ')' + text[(tmp.span()[1]+1):-1] + '\n'
                text = text.replace("_", "0")       # rmarkdown do not support _ in labels
                
            # equation ref
            tmp = re.search('\\\\eqref{[\S\s]*?}', text)
            if not tmp is None:
                num_ref = len(re.findall('\\\\eqref{[\S\s]*?}', text))
                ii = 0
                while ii < num_ref:
                    tmp_end = re.search('\\\\eqref{[\S\s]*?}', text)
                    tmp_begin = re.search('\\\\eqref{', text)
                    tmp_cont = text[(tmp_begin.span()[0]+7): (tmp_end.span()[1]-1)]
                    tmp_cont = tmp_cont.replace("_", "0")
                    text = text[0:(tmp_begin.span()[0]-1)] + '\@ref(' + tmp_cont + ')' + text[(tmp_end.span()[1]+1):-1]
                    # print(text)
                    ii += 1
                
            file_out.write(text)

    file_in.close()
    file_out.close()

    # os.remove(file_path + file_name)
    # os.rename(file_path + 'tmp.html', file_path + file_name)
