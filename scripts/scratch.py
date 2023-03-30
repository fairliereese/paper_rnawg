import glob

for f in glob.glob('*sam'):
    ifile = open(f, 'r')
    ofile = open(f+'_new', 'w')

    new_chroms = []

    for line in ifile:
        if line.startswith('@'):
            ofile.write(line)
        else:
            line = line.strip().split('\t')
            line[2] = line[2].split('_', maxsplit=1)[-1]
            new_chroms.append(line[2])
            line = '\t'.join(line)+'\n'
            ofile.write(line)
