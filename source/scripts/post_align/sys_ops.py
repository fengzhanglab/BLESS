import csv
import sys
import time

class ErrorLog:
    def __init__(self):
        self.my_times = []
        self.my_exceptions = []
    def append(self, time_exception):
        #appends list of lists with time and exception as first and second elements, respectively
        self.my_times = [time_exception[i][0] for i in range(len(time_exception))]
        self.my_exceptions = [time_exception[i][1] for i in range(len(time_exception))]
        
def exitProgram():
    print "Press any key to exit ..."
    raw_input()
    sys.exit()
            
def throw_exception(this_input):
    #throws exception this_input[0] to file-name this_input[1], if this_input[1] exists, or errorlog.csv otherwise

    if(type(this_input)==list and len(this_input)==2):
        errorphrase = this_input[0]
        errorlog_filename = this_input[1]
    else:
        if(type(this_input)==list):
            errorphrase = this_input[0]
        else:
            errorphrase = this_input
        errorlog_filename = "errorlog.csv"        

    my_datetime = time.strftime("%Y/%m/%d %H:%M:%S")

    print errorphrase
    print my_datetime
    
    my_errorlog = ErrorLog()

    try:
        with open(errorlog_filename,'r') as csvfile:
            errorlog_reader = csv.reader(csvfile, delimiter='|')
            print "Opened " + errorlog_filename
            my_errorlog.append([row for row in errorlog_reader])
    except:
        print "Creating new error-log file, " + errorlog_filename + " as none exists yet ..."

    my_errorlog.my_times.append(my_datetime)
    my_errorlog.my_exceptions.append(" " + errorphrase)

    with open(errorlog_filename,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter='|')
        mywriter.writerows([[my_errorlog.my_times[i], my_errorlog.my_exceptions[i]] for i in range(len(my_errorlog.my_times))])

    csvfile.close()
    print my_datetime + ", " + errorphrase

def throw_status(this_input):
    #throws status this_input[0] to file-name this_input[1], if this_input[1] exists, or statuslog.csv otherwise

    if(type(this_input)==list and len(this_input)==2):
        statusphrase = this_input[0]
        statuslog_filename = this_input[1]
    else:
        if(type(this_input)==list):
            statusphrase = this_input[0]
        else:
            statusphrase = this_input
        statuslog_filename = "statuslog.csv"        

    my_datetime = time.strftime("%Y/%m/%d %H:%M:%S")

    print statusphrase
    print my_datetime
    
    my_statuslog = ErrorLog()

    try:
        with open(statuslog_filename,'r') as csvfile:
            statuslog_reader = csv.reader(csvfile, delimiter='|')
            print "Opened " + statuslog_filename
            my_statuslog.append([row for row in statuslog_reader])
    except:
        print "Creating new status-log file, " + statuslog_filename + " as none exists yet ..."

    my_statuslog.my_times.append(my_datetime)
    my_statuslog.my_exceptions.append(" " + statusphrase)

    with open(statuslog_filename,'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter='|')
        mywriter.writerows([[my_statuslog.my_times[i], my_statuslog.my_exceptions[i]] for i in range(len(my_statuslog.my_times))])

    csvfile.close()
    print my_datetime + ", " + statusphrase
