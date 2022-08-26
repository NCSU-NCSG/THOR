import sys

comp_list_file = open("reg_tests_comp_list.txt","r")
for lline in comp_list_file:
    comp_files=lline.split()
    f1 = open(comp_files[0])
    f2 = open(comp_files[1])
    print(comp_files[0])
    i = 0
    for line1 in f1:
        i += 1
        for line2 in f2:
            # matching line1 from both files
            if line1 != line2:
                line1=line1.replace(',',' ')
                line1=line1.replace("'",' ')
                line2=line2.replace(',',' ')
                line2=line2.replace("'",' ')
                line1split=line1.split(' ')
                line2split=line2.split(' ')
                for ii in range(len(line1split)):
                    if line1split[ii] != line2split[ii]:
                        if abs(float(line1split[ii])-float(line2split[ii])) >= 1.0e-10:
                            print("-----------------------------------DIFFERENCE FOUND!-----------------------------------")
                            print("\tFile 1:", line1, end='')
                            print("\tFile 2:", line2, end='')
            break
    f1.close()
    f2.close()
