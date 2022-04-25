import glob, os

files = []
for filename in glob.iglob('./refractiveindex.info-database/database/data/**', recursive=True):
    if os.path.isfile(filename): # filter dirs
        with open(filename) as f:
            contents = f.read()
            # if contents.count('- type:')>1:
            #     files.append(filename.replace('./refractiveindex.info-database/database/data/',''))
            if 'tabulated n' in contents:
                files.append(filename.replace('./refractiveindex.info-database/database/data/',''))

print(files)
print(len(files))
with open('tabulated.txt', 'w') as f:
    f.write(str(files))
