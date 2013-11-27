import subprocess
p = subprocess.Popen('python child.py', shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
stdout, stderr = p.communicate()

print stdout == ""
print stderr == ""

#print stdout
#print stderr
