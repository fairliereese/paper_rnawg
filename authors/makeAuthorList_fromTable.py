#!/usr/bin/python
# <Script name>
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse
import os
import pdb
import csv
import io
import codecs

#############
# CONSTANTS #
#############

#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print("%s option not supplied" % option)
            self.print_help()
            sys.exit(1)

# from https://stackoverflow.com/questions/5004687/python-csv-dictreader-with-utf-8-data
#class UnicodeDictReader(csv.DictReader, object):
#    def next(self):
#        try:
#            row = super(UnicodeDictReader, self).next()
#            return {unicode(key, 'utf-8'): unicode(value, 'utf-8') for key, value in row.iteritems()}
#	except:
#            pdb.set_trace()

###############
# END CLASSES #
###############

########
# MAIN #
########
def main():

    opt_parser = OptionParser()

    # Add Options. Required options should have default=None
    opt_parser.add_option("-i",
                          dest="input",
                          type="string",
                          help="""Input table with columns "First, Middle, Last,
                                  City, State, Country, Institution."
                                  Tab-delmited. Order of author list must be
                                  order in table. Table should be pretty clean
                                  and authors with multiple affiliations are
                                  listed consecutively with each affiliation""",
                         default=None)
    opt_parser.add_option("-d",
                          dest="use_department",
                          action="store_true",
                          help="""Input table also has department field;
                                  therefore, also include department""",
                         default=False)
    opt_parser.add_option("--rtf",
                          dest="rtf",
                          type="string",
                          help="""Output file. Formatted authors and affiliations in RTF
                                  format to include superscripts""",
                         default=None)
    opt_parser.add_option("--no_affiliations",
                          dest="no_affil",
                          action="store_true",
                          help="""Only print out names in order. Do not print
                                  out affiliations""",
                         default=False)

    (options, args) = opt_parser.parse_args()

    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("--rtf")

    use_department = options.use_department

    name_dict = {}
    num2affil = {}
    affil2num = {}
    affil_ctr = 0

    name_order = []

#   infile = codecs.open(options.input, mode='r', encoding='utf-8')


#    csv_reader = UnicodeDictReader(infile, '\t'iter="\t")

    infile = open(options.input, mode='r')

#    out_rtf = open(options.rtf, "w")
    out_rtf = codecs.open(options.rtf, "w", encoding="utf-8")

    csv_reader = UnicodeDictReader(infile, delimiter=",")

    no_affil = options.no_affil

    for row in csv_reader:

        name, affil = parseLine(row, use_department)

        if affil in affil2num:
            affil_num = affil2num[affil]
        else:
            affil_ctr += 1
            num2affil[affil_ctr] = affil
            affil2num[affil] = affil_ctr
            affil_num = affil_ctr

        if len(name_order) == 0:
            name_order.append(name)
            name_dict[name] = [affil_num]
        elif name == name_order[-1]:
            name_dict[name].append(affil_num)
        else:
            name_order.append(name)
            name_dict[name] = [affil_num]

    infile.close()


    # print to rtf
    out_rtf.write("{\\rtf\\ansi\\ansicpg1252\\deff0 {\\fonttbl {\\f0 Arial;}}\\fs24\n")

    #print author list
    out_rtf.write("{\\pard\\sa180\n")

    name_str = ""
    for name in name_order:
        name_dict[name] = sorted(name_dict[name])
        if no_affil:
            name_str += ", %s" % name
        else:
            name_str += ", %s{\\super %s}" % (name,
                                          ",".join(map(repr,name_dict[name])))

    name_str = name_str.lstrip(", ")

    # end author block
    out_rtf.write("%s\n\\par}\n" % name_str)

    # Start affil block
    out_rtf.write("{\\pard\\sa180\n")

    affil_str = ""
    affil_list = num2affil.keys()
    affil_list = sorted(affil_list)

    # Remove redundancy of institute information
    if use_department:
        num2institute_flag = {} # helps to flag when to print(the institute
                                # along with department
        # With departments, the affiliation is a tuple
        for num in affil_list:
            num2institute_flag[num] = True
            if num == 1:
                continue
            if num2affil[num][1] == num2affil[num-1][1]:
                num2institute_flag[num-1] = False

    for num in affil_list:
        if use_department:
            if num2institute_flag[num]:
                # Checking for blank departments
                if num2affil[num][0] == "":
                    # If no department informaiton, should at least have
                    # institutional information
                    this_affil = num2affil[num][1]
                else:
                    this_affil = ", ".join(num2affil[num])
            else:
                # Just use department flag
                this_affil = num2affil[num][0]
        else:
            this_affil = num2affil[num]

        affil_str += ", {\\super %d}%s" % (num, this_affil)

    affil_str = affil_str.lstrip(", ")

    out_rtf.write("%s\n\\par}\n}" % affil_str)

    out_rtf.close()



    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def parseLine(row, use_department):

    m = row['Middle']
    if m == "":
        middle = " "
    else:
        if len(m)==1:
            middle = ' {}. '.format(m)
        else:
            middle = ' {} '.format(m)
    if len(row['First']) == 1:
        first = '{}.'.format(row['First'])
    else:
        first = row['First']

    name = '{}{}{}'.format(first, middle, row['Last'])
    # import pdb
    # pdb.set_trace()
    # if row["Middle"] == "":
    #     name = "%s %s" % (row["First"], row["Last"])
    # else:
    #     name = "%s %s. %s" % (row["First"], row["Middle"], row["Last"])

    affil_list = []
    if row["Institution"] != "":
        affil_list.append(row["Institution"])
    if row["City"] != "":
        affil_list.append(row["City"])

    try:
        if row["Zipcode"] != "":
            affil_list.append(row["Zipcode"])
    except:
        print("Warning: No Zipcode column")

    if row["Country"] != "":
        affil_list.append(row["Country"])

    if use_department:
        try:
            department = row["Department"]
        except:
            print("Can't find \'Department\' column.")
            sys.exit(1)

        affil = (department, ", ".join(affil_list))

    else:
        affil = ", ".join(affil_list)

#   affil = "%s, %s, %s, %s, %s" % (row["Institution"],
#                               row["City"],
#                               row["State"],
#                               row["Zipcode"],
#                               row["Country"])

    return name, affil

def UnicodeDictReader(utf8_data, **kwargs):
    csv_reader = csv.DictReader(utf8_data, **kwargs)
    for row in csv_reader:
# 	pdb.set_trace()
#        yield {codecs.decode(key, 'utf-8'):codecs.decode(value,'utf-8') for key, value in row.iteritems()}
#        pdb.set_trace()
        yield {key:value for key, value in row.items()}



#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()
