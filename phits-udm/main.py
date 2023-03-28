import os, sys, glob
import diff


# ======================================================================
SupportedVersions=["3.30"]


# ======================================================================
DIR_ori="../src"      # Original directory
DIR_udm="../src-udm"  # New directory for phits-udm


# ======================================================================
if os.path.exists(DIR_udm):
    sys.exit("@@@ ERROR: 'src-udm' already exists.")


# ======================================================================
# get version
with open(DIR_ori+"/main.f") as f:
    lines=f.read().split('\n')

for l in lines:
    if "Version" in l:
        v=l.split()[4]
        break

if v not in SupportedVersions:
    print("Version error")
    print("This python script does not support your version, "+v+".")
    print("Supported versions are:")
    for sv in SupportedVersions:
        print("  "+sv)
    sys.exit()


if v=="3.30":
    fL=diff.v330.filename
    bL=diff.v330.before
    aL=diff.v330.after


# ======================================================================
os.system("cp -r {} {}".format(DIR_ori,DIR_udm))
os.system("cp udm_Parameter.f {}/.".format(DIR_udm))
# os.system("cp sample-code-1/udm_Manager.f90 {}/.".format(DIR_udm))
# os.system("cp sample-code-1/udm_int_sample*.f90  {}/.".format(DIR_udm))
# os.system("cp sample-code-1/udm_part_sample*.f90 {}/.".format(DIR_udm))
os.chdir(DIR_udm)
os.system("make clean")


# ======================================================================
# replace
def check(f,b,a):
    if not os.path.exists(f):
        sys.exit("not found: "+f)
    with open(f) as FILE:
        s=FILE.read()
    if s.count(b)!=1:
        print("@@@ N =",s.count(b))
        print(f)
        print(b)
        print(a)
        sys.exit("ERROR")

for i,[f,b,a] in enumerate(zip(fL,bL,aL)):
    check(f,b,a)
    with open(f    ) as FILE: s=FILE.read()
    s=s.replace(b,a)
    with open(f,"w") as FILE: FILE.write(s)
    print(i+1,f)


print("!! Completed !!")





