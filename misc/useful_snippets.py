
##############################
# Insert code from funcdefs.py into a script 
# >Gives, e.g., functions defined in funcdefs.py access to global-scope variables in script.py
insertme = './testjerk.py' 
print(f"Inserting {insertme:s} into this script ...")
exec(open(insertme).read())
