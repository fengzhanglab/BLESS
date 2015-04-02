import sys
import csv
import sys_ops

def generate_sample_jobfile(param_filename,process_filename,jobfile_filename):
    
    [f_barcodes, r_barcodes, groups] = get_params(param_filename)
    
    process_list_file = []
    try:
        with open(process_filename, 'rb') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            process_file = [row for row in csvreader]
    except:
        sys_ops.throw_exception("Could not open process-list " + process_filename)
        sys.exit(1)
    
    job_list = []
    for row in process_file:
        print row
        #columns in process_list_file must be, in order:
        f_barcode_checkList = [f_barcode in str(row[0]) for f_barcode in f_barcodes]
        r_barcode_checkList = [r_barcode in str(row[1]) for r_barcode in r_barcodes]
        #f_barcode_checkList2 = [f_barcode in str(row[3]) for f_barcode in f_barcodes]
        #r_barcode_checkList2 = [r_barcode in str(row[4]) for r_barcode in r_barcodes]

        #NOTE: in each of the above cases, the reader assumes the stored sequences in the function get_sequences() are necessary in their entirety to exist as substrings of the sequences in the process-list
        #if((True in f_barcode_checkList) and (True in r_barcode_checkList) and (True in f_barcode_checkList2) and (True in r_barcode_checkList2)):
        if 0:
            print row

            if((len(row[0])==0) and (len(row[1])==0)):             
                job_list.append([str(0), str(-1)+'_'+str(-1)]) 
                job_list.append([str(0), str(-1)+'_'+str(-1)])                        
            elif(len(row[0])==0):             
                job_list.append([str(0), str(-1)+'_'+str(r_barcode_checkList.index(True))]) 
                job_list.append([str(0), str(-1)+'_'+str(r_barcode_checkList2.index(True))])
            elif(len(row[1])==0): 
                job_list.append([str(0), str(f_barcode_checkList.index(True))+'_'+str(-1)])
                job_list.append([str(0), str(f_barcode_checkList.index2(True))+'_'+str(-1)])                                       
            else:
                job_list.append([str(0), str(f_barcode_checkList.index(True))+'_'+str(r_barcode_checkList.index(True))])
                job_list.append([str(0), str(f_barcode_checkList2.index(True))+'_'+str(r_barcode_checkList2.index(True))])
        elif ((True in f_barcode_checkList) and (True in r_barcode_checkList)):
            print row

            if((len(row[0])==0) and (len(row[1])==0)):             
                job_list.append([str(0), str(-1)+'_'+str(-1)])                        
            elif(len(row[0])==0):             
                job_list.append([str(0), str(-1)+'_'+str(r_barcode_checkList.index(True))])
            elif(len(row[1])==0): 
                job_list.append([str(0), str(f_barcode_checkList.index(True))+'_'+str(-1)])                                       
            else:
                job_list.append([str(0), str(f_barcode_checkList.index(True))+'_'+str(r_barcode_checkList.index(True))])            
                #include 0 as first element to indicate that the currently-written job has NOT been completed
        elif((True in f_barcode_checkList) or (True in r_barcode_checkList) or (True in f_barcode_checkList2) or (True in r_barcode_checkList2)):
            sys_ops.throw_exception("Process list row " + str(row) + " contains mixture of identifiable and unidentifiable sequences. Exiting.")
            sys.exit(1)

        print job_list            
    with open(jobfile_filename,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        for job in job_list:
            mywriter.writerow(job)
            
    return 

def generate_replicate_jobfile(param_filename,process_filename,jobfile_filename):
    
    [f_barcodes, r_barcodes, groups] = get_params(param_filename)
    
    process_list_file = []
    try:
        with open(process_filename, 'rb') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            process_file = [row for row in csvreader]
    except:
        sys_ops.throw_exception("Could not open process-list " + process_filename)
        sys.exit(1)
    
    groupids = []
    job_list = {}
    for row in process_file:
        print row
        #columns in process_list_file must be, in order:
        f_barcode_checkList = [f_barcode in str(row[0]) for f_barcode in f_barcodes]
        r_barcode_checkList = [r_barcode in str(row[1]) for r_barcode in r_barcodes]
        group_checkList = [group in str(row[2]) for group in groups]

        #NOTE: in each of the above cases, the reader assumes the stored sequences in the function get_sequences() are necessary in their entirety to exist as substrings of the sequences in the process-list
        if((True in f_barcode_checkList) and (True in r_barcode_checkList) and (True in group_checkList)):
            print row
            if (hash(str(row[2])) not in job_list):
                groupids.append(hash(str(row[2])))
                job_list[hash(str(row[2]))] = []
                job_list[hash(str(row[2]))].append(str(0))

            if((len(row[0])==0) and (len(row[1])==0)):             
                job_list[hash(str(row[2]))].append(str(-1)+'_'+str(-1))                        
            elif(len(row[0])==0):             
                job_list[hash(str(row[2]))].append(str(-1)+'_'+str(r_barcode_checkList.index(True)))
            elif(len(row[1])==0): 
                job_list[hash(str(row[2]))].append(str(f_barcode_checkList.index(True))+'_'+str(-1))                                       
            else:
                job_list[hash(str(row[2]))].append(str(f_barcode_checkList.index(True))+'_'+str(r_barcode_checkList.index(True)))

                #include 0 as first element to indicate that the currently-written job has NOT been completed
        elif((True in f_barcode_checkList) or (True in r_barcode_checkList) or (True in group_checkList)):
            sys_ops.throw_exception("Process list row " + str(row) + " contains mixture of identifiable and unidentifiable sequences. Exiting.")
            sys.exit(1)

        print job_list            
    with open(jobfile_filename,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        for groupid in groupids:
            print groupid
            mywriter.writerow(job_list[groupid])
            
    return

def generate_junction_jobfile(param_filename,process_filename,jobfile_filename):
    
    [f_barcodes, r_barcodes, groups] = get_params(param_filename)
    
    process_list_file = []
    try:
        with open(process_filename, 'rb') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            process_file = [row for row in csvreader]
    except:
        sys_ops.throw_exception("Could not open process-list " + process_filename)
        sys.exit(1)
    
    job_list = []
    for row in process_file:
        print row
        #columns in process_list_file must be, in order:
        f_barcode_checkList = [f_barcode in str(row[0]) for f_barcode in f_barcodes]
        r_barcode_checkList = [r_barcode in str(row[1]) for r_barcode in r_barcodes]
        f_barcode_checkList2 = [f_barcode in str(row[3]) for f_barcode in f_barcodes]
        r_barcode_checkList2 = [r_barcode in str(row[4]) for r_barcode in r_barcodes]

        #NOTE: in each of the above cases, the reader assumes the stored sequences in the function get_sequences() are necessary in their entirety to exist as substrings of the sequences in the process-list
        if((True in f_barcode_checkList) and (True in r_barcode_checkList) and (True in f_barcode_checkList2) and (True in r_barcode_checkList2)):
            print row

            if((len(row[0])==0) and (len(row[1])==0)):             
                job_list.append([str(0), str(-1)+'_'+str(-1), str(-1)+'_'+str(-1)])                        
            elif(len(row[0])==0):             
                job_list.append([str(0), str(-1)+'_'+str(r_barcode_checkList.index(True)), str(-1)+'_'+str(r_barcode_checkList2.index(True))])
            elif(len(row[1])==0): 
                job_list.append([str(0), str(f_barcode_checkList.index(True))+'_'+str(-1), str(f_barcode_checkList.index2(True))+'_'+str(-1)])                                       
            else:
                job_list.append([str(0), str(f_barcode_checkList.index(True))+'_'+str(r_barcode_checkList.index(True)), str(f_barcode_checkList2.index(True))+'_'+str(r_barcode_checkList2.index(True))])

                #include 0 as first element to indicate that the currently-written job has NOT been completed
        elif((True in f_barcode_checkList) or (True in r_barcode_checkList) or (True in f_barcode_checkList2) or (True in r_barcode_checkList2)):
            sys_ops.throw_exception("Process list row " + str(row) + " contains mixture of identifiable and unidentifiable sequences. Exiting.")
            sys.exit(1)

        print job_list            
    with open(jobfile_filename,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        for job in job_list:
            mywriter.writerow(job)
            
    return    

def get_params(readparam_fileName):
    readparam_file = []
    try:
        with open(readparam_fileName, 'rb') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            readparam_file = [row for row in csvreader]
    except:
        sys_ops.throw_exception("Could not open read-params " + readparam_fileName)
        sys.exit(1)
##############################

    f_barcodes = readparam_file[0]
    if '' in f_barcodes:
        temp_index = f_barcodes.index('')
        f_barcodes = f_barcodes[0:temp_index]
        for i in range(0,len(f_barcodes)):
            f_barcodes[i] = str(f_barcodes[i])    

    r_barcodes = readparam_file[1]
    if '' in r_barcodes:
        temp_index = r_barcodes.index('')
        r_barcodes = r_barcodes[0:temp_index]
        for i in range(0,len(r_barcodes)):  
            r_barcodes[i] = str(r_barcodes[i]) 

    groups = readparam_file[2]
    if '' in groups:
        temp_index = groups.index('')
        groups = groups[0:temp_index] 
        for i in range(0,len(groups)):
            groups[i] = str(groups[i])      

##############################
    return [f_barcodes, r_barcodes, groups]

def get_next_job(jobfile_filename,fromVal,toVal): #opens job_file, returns barcode/primer pair for first row containing fromVal in first column
    f_barcode_ind = -1
    r_barcode_ind = -1
    new_job_file = []
    with open(jobfile_filename, 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        new_job_file =  [row for row in csvreader]

    job = []
    with open(jobfile_filename,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',')
        found = False
        for row in new_job_file:
            if(str(row[0])==str(fromVal) and (not found)):
                row[0] = str(toVal)
                job = row[1:len(row)]
                found = True
            mywriter.writerow(row)
    
    return job

