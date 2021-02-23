#This python script was created on a Linux Ubuntu 18.04 operating syste using PyCharm community edition 
#Parse from phobius will allow parsing the sequences obtained by the analysis performed with phobius. 
#As it is now, the script will allow you to parse all the receptors that contain from 4 to 9 transmembrane domains. This value can be changed 
#by modifying the range in line 24.
#The input needed for this script is a CSV file produced by phobius. 
#The output will be a list containing only the ID of the sequences that have between 4 to 9 transmembrane domains. 


import csv

# creating the list1
tm_morethan4tm = []
with open('phobiusfortestafterdeletedrubbish.csv', 'r') as project_csv:
    # Delimiter can be changed according to what it is, in this case i have whitesoace
    # Which is the default, phobius will produce a whitespace delimited file, you can change the whitespace for
    # white space in this part
    csv_dict_reader = csv.DictReader(project_csv, delimiter=' ')

    for row in csv_dict_reader:
        # the number of transmembrane domains should be above certain threshold
        tm = int(row['TM'])
        # here you can add any if condition in this case I'm searching for all the sequences that show
        # more than 4 transmembrane domains
        if tm in range(4, 9):
            # write a new CSV appending the data that is coming from that is coming back from the CSV.
            tm_morethan4tm.append({
                'TM': row['TM'],
                'ID': row['SEQUENCE']
            })
# This line is to create a new file and help to 'w' Write on it. Newline is needed because if not pythonwill leave
# a blank line after each row
# the as is
with open('tmfinal.csv', 'w', newline='') as tm_morethan4tms:
    # Fieldnames can be changed to put the ID or TM first
    fieldnames = ['ID', 'TM']
    csv_dict_writer = csv.DictWriter(tm_morethan4tms, fieldnames=fieldnames)
    csv_dict_writer.writeheader()

    for tm in tm_morethan4tm:
        csv_dict_writer.writerow(tm)
