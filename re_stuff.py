import re

string = "11:72297209_C/T_exm938858"

location = re.search('\S+:\S+_[ACTG]/[ACTG]', string)

print location.group()