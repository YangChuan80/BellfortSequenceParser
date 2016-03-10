# Bellfort Sequence Parser

## Modules

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk
import tkinter.font as tkf
from tkinter import messagebox
from tkinter import filedialog
import threading
import time
import os
import shutil

## Helper Functions
### Reverse Complement

def reverseComplement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    rc_sequence=''
    for s in sequence:
        rc_sequence = complement[s] + rc_sequence
    return rc_sequence

### FASTQ File Browse

def buttonBrowseFASTQ():
    global filenameFASTQ, indicator_preprocess
    
    try:
        filenameFASTQ = filedialog.askopenfilename(filetypes=(('FASTQ files', '*.fastq'), 
                                                              ('All files', '*.*')))
        text_fileFASTQ.delete('1.0', tk.END)
        text_fileFASTQ.insert('1.0', filenameFASTQ.split('/')[-1])
        
        # Reset the progress bar///////////////
        progressbar['value'] = 0
        progressbar_loadFASTQ['value'] = 0
        
        # Reset the percentage
        text_percentage.delete('1.0', tk.END)
        text_percentage.insert('1.0', str('0%'))
        
        indicator_preprocess = 0
    except:
        filenameFASTQ = ''   

### FASTQ File Load

def loadFASTQ():
    global reads
    
    start_time = time.time()    
    
    f = open(filenameFASTQ)

    reads = []

    try:
        while 1:
            name = f.readline().rstrip()
            sequence = f.readline().rstrip()
            f.readline()
            quality = f.readline().rstrip()

            if len(name) == 0:
                break

            union = name, sequence

            reads.append(union)           

        end_time = time.time()
        delta_time = end_time - start_time

        text_time.delete('1.0', tk.END)
        text_time.insert('1.0', str(delta_time))  

        text_readNum.delete('1.0', tk.END)
        text_readNum.insert('1.0', str(len(reads)))  

    except:
        messagebox.showwarning("File Loading Failed", 
                               "Sorry, file loading failed! Please check the file format.")
    f.close()

def start_loadFASTQ_thread(event):
    global loadFASTQ_thread
    
    if filenameFASTQ != '':
        loadFASTQ_thread = threading.Thread(target=loadFASTQ)
        loadFASTQ_thread.daemon = True

        progressbar_loadFASTQ.start(10)
        loadFASTQ_thread.start()
        root.after(20, check_loadFASTQ_thread)
    else:
        messagebox.showwarning("No File", 
                               "Sorry, no file loaded! Please choose FASTQ file first.")

def check_loadFASTQ_thread():
    if loadFASTQ_thread.is_alive():
        progressbar_loadFASTQ.start(10)
        root.after(20, check_loadFASTQ_thread)
    else:
        progressbar_loadFASTQ.stop()
        progressbar_loadFASTQ['value']=100
        messagebox.showinfo("FASTQ File Loaded", "FASTQ file successfully loaded!")

### Divide FASTQ File

def divideFASTQ():    
    start_time = time.time() 

    gotten = text_readNumDivided.get('1.0', tk.END)
    readNumDivided = int(gotten.rstrip())

    if os.path.exists(filenameFASTQ+'.folder'):           
        # Remove the folder previously made:
        shutil.rmtree(filenameFASTQ+'.folder')

    # Make a new one:
    os.makedirs(filenameFASTQ+'.folder')    

    line_num = 0
    file_no = 1

    f_input = open(filenameFASTQ)        
    f_output = open(filenameFASTQ+'.folder/' + 'Reads_Slice_No_' + str(file_no) + '.fastq', 'w') 

    while 1:
        # Input ///////////////////////////////////
        name = f_input.readline()
        sequence = f_input.readline()
        f_input.readline()
        quality = f_input.readline()

        if len(name) == 0:
            break            

        # Output ////////////////////////////////////

        f_output.write(name)
        f_output.write(sequence)
        f_output.write('+\n')
        f_output.write(quality)      

        line_num += 1

        if line_num == readNumDivided:                
            f_output.close() 
            file_no += 1
            f_output = open(filenameFASTQ+'.folder/' + 'Reads_Slice_No_' + str(file_no) + '.fastq', 'w')
            line_num = 0                

    end_time = time.time()
    delta_time = end_time - start_time

    text_time.delete('1.0', tk.END)
    text_time.insert('1.0', str(delta_time))  

    f_input.close()
    f_output.close()

def start_divideFASTQ_thread(event):
    global divideFASTQ_thread
    
    if filenameFASTQ != '':
        divideFASTQ_thread = threading.Thread(target=divideFASTQ)
        divideFASTQ_thread.daemon = True

        progressbar_loadFASTQ.start(10)
        divideFASTQ_thread.start()
        root.after(20, check_divideFASTQ_thread)
    else:
        messagebox.showwarning("No File", 
                               "Sorry, no file loaded! Please choose FASTQ file first.")

def check_divideFASTQ_thread():
    if divideFASTQ_thread.is_alive():
        progressbar_loadFASTQ.start(10)
        root.after(20, check_divideFASTQ_thread)
    else:
        progressbar_loadFASTQ.stop()
        progressbar_loadFASTQ['value']=100
        messagebox.showinfo("FASTQ File Divided", "FASTQ file has been successfully divided!")

### Preprocess

def preprocessFASTQ():
    global reads, indicator_preprocess, kmer_dict_reads
    
    try:
        num = len(reads)   
        indicator_preprocess = 0
        gain = 50/num

        gotten = text_sequence_len.get('1.0', tk.END)
        k = int(gotten.rstrip())
        
        if k > len(reads[0][1]):
            messagebox.showwarning("Target Sequence Length Error", 
                                   "Sorry, the target sequence length is more than read length. Please check.")
        elif k < 3:
            messagebox.showwarning("Sequence Too Short", 
                                   "Sorry, the target sequence length is too short which will make the program running slowly. Please check.")
        elif filenameSequences == '':
            messagebox.showwarning("No Sequences Loaded", 
                                   "Sorry, no sequences loaded! Please load sequences first.")
        else:
            kmer_dict_reads = {}

            start_time = time.time()

            for read in reads:
                for i in range(len(read[1])-k+1):
                    kmer_dict_reads[read[1][i:i+k]] = set()
                indicator_preprocess += gain 

            for read in reads:
                for i in range(len(read[1])-k+1):
                    kmer_dict_reads[read[1][i:i+k]].add(read)
                indicator_preprocess += gain
                
            indicator_progress = 100
            
            # Add MatchAll Here ///////////////////////////////////////////////////////
            matchAll()

            end_time = time.time()
            delta_time = end_time - start_time

            text_time.delete('1.0', tk.END)
            text_time.insert('1.0', str(delta_time)) 
            
            messagebox.showinfo("Preprocess FASTQ & Count Matched Sequences Completed", 
                                "Current FASTQ preprocess & matched sequence counts successfully completed!")

    except NameError:
        messagebox.showwarning("No FASTQ File Loaded", 
                               "Sorry, no loaded FASTQ file found! Please load FASTQ file first.")

def start_preprocess_thread(event):
    global preprocess_thread, indicator_preprocess
    preprocess_thread = threading.Thread(target=preprocessFASTQ)
    preprocess_thread.daemon = True
    
    progressbar['value'] = indicator_preprocess
    text_percentage.delete('1.0', tk.END)
    text_percentage.insert('1.0', str(int(indicator_preprocess))+'%')
    
    preprocess_thread.start()
    root.after(20, check_preprocess_thread)

def check_preprocess_thread():
    if preprocess_thread.is_alive():
        progressbar['value'] = indicator_preprocess
        text_percentage.delete('1.0', tk.END)
        text_percentage.insert('1.0', str(int(indicator_preprocess))+'%')
        
        root.after(20, check_preprocess_thread)

### Match All

def matchAll():
    global  kmer_dict_reads, indicator_matchAll, df
    
    try:
        len(kmer_dict_reads)    
        num = len(df)
        
        if num == 0:
            messagebox.showwarning("No Sequences Loaded", 
                                   "Sorry, no sequences loaded! Please load sequences first.")
        else:    
            indicator_matchAll = 0
            gain = 1000000/num

            start_time = time.time()

            arr = np.array(df)

            for i in range(len(arr)):
                key1 = arr[i,2]
                key2 = reverseComplement(key1)
                
                try:
                    n1 = len(kmer_dict_reads[key1])
                except KeyError:
                    n1 = 0
                    
                try:
                    n2 = len(kmer_dict_reads[key2])
                except KeyError:
                    n2 = 0
                    
                arr[i, 4] += n1 + n2
                arr[i, 5] += 1
                
                indicator_matchAll += gain

            df = pd.DataFrame(arr, columns = ['gene_id', 'UID', 'seq', 'Reserved', 'Count', 'Tag'])
            #df = df.set_index('UID', drop=False) 

            end_time = time.time()
            delta_time = end_time - start_time

            text_time.delete('1.0', tk.END)
            text_time.insert('1.0', str(delta_time))           

    except NameError:
        messagebox.showwarning("No FASTQ Preprocessed or No Sequences Loaded", 
                               "Sorry, no FASTQ preprocess implemented or no sequences file loaded! Please preprocess FASTQ or load sequences first.")    

def start_matchAll_thread(event):
    global matchAll_thread, indicator_matchAll
    matchAll_thread = threading.Thread(target=matchAll)
    matchAll_thread.daemon = True
    
    progressbar['value'] = indicator_matchAll
    
    matchAll_thread.start()
    root.after(20, check_matchAll_thread)

def check_matchAll_thread():
    if matchAll_thread.is_alive():
        progressbar['value'] = indicator_matchAll
        
        root.after(20, check_matchAll_thread)
    else:
        messagebox.showinfo("Matching Completed", 
                                "Counting of sequences matched successfully completed!")

### Match Single

def buttonMatch():
    gotten = text_sequence.get('1.0', tk.END)
    p1 = gotten.rstrip()    
    p2 = reverseComplement(p1)
    
    if p1 == '' or p2 == '':
        messagebox.showwarning("No Sequence Found", 
                               "Sorry, no sequence found in the text blank above! Please check the sequence.")
    else:
        try:
            len(kmer_dict_reads)
            try:
                n1 = len(kmer_dict_reads[p1])
            except KeyError:
                n1 = 0
            
            try:
                n2 = len(kmer_dict_reads[p2])
            except KeyError:
                n2 = 0
                
            count = n1 + n2
                
            text_count.delete('1.0', tk.END)
            text_count.insert('1.0', str(count))
            
        except NameError:
            messagebox.showwarning("No FASTQ Preprocessed", 
                                   "Sorry, no FASTQ preprocess implemented! Please preprocess FASTQ first.")

### File of Target Sequence Load

def loadSequences():
    global filenameSequences, df, recordNum
    
    progressbar_loadSequences['value'] = 0
    try:
        filenameSequences = filedialog.askopenfilename(filetypes=(('Comma-Separated (CSV) text file', '*.csv'), ('All files', '*.*')))
        text_fileSequences.delete('1.0', tk.END)
        text_fileSequences.insert('1.0', filenameSequences.split('/')[-1])
    except:
        filenameSequences = ''    
   
    if filenameSequences == '':
        messagebox.showwarning("No File", "Sorry, no file chosen! Please choose file of sequences first.")
    else:        
        try:
            start_time = time.time()
            
            df = pd.read_csv(filenameSequences)
            df['count'] = 0
            df['tag'] = 0
            #df = df.set_index('UID', drop=False)  
            
            recordNum = len(df)
            
            progressbar_loadSequences['value'] = 100
            
            end_time = time.time()
            delta_time = end_time - start_time
                       
            text_time.delete('1.0', tk.END)
            text_time.insert('1.0', str(delta_time))
            
            text_recordNum.delete('1.0', tk.END)
            text_recordNum.insert('1.0', str(recordNum))
            
            messagebox.showinfo("File of Sequences Loaded", "File of sequences successfully loaded!")        
        except:
            messagebox.showwarning("File Loading Failed", "Sorry, file loading failed! Please check the file format.")    

### Load Half-matched Sequences

def buttonLoadHalfMatchedSequences():
    global df, recordNum, filenameSequences
    
    progressbar_loadSequences['value'] = 0
    try:
        filenameSequences = filedialog.askopenfilename(filetypes=(('Comma-Separated (CSV) text file', '*.csv'), ('All files', '*.*')))
        text_fileSequences.delete('1.0', tk.END)
        text_fileSequences.insert('1.0', filenameSequences.split('/')[-1])
    except:
        filenameSequences = ''    
        
    if filenameSequences == '':
        messagebox.showwarning("No File", "Sorry, no file chosen! Please choose file of sequences first.")
    else:        
        try:
            start_time = time.time()
            
            df = pd.read_csv(filenameSequences)   
            df = df.set_index('Unnamed: 0', drop=True)  
            
            recordNum = len(df)
            
            progressbar_loadSequences['value'] = 100
            
            end_time = time.time()
            delta_time = end_time - start_time
                       
            text_time.delete('1.0', tk.END)
            text_time.insert('1.0', str(delta_time))
            
            text_recordNum.delete('1.0', tk.END)
            text_recordNum.insert('1.0', str(recordNum))
            
            messagebox.showinfo("File of Half Matched Sequences Loaded", "File of half matched sequences successfully loaded!")        
        except:
            messagebox.showwarning("File Loading Failed", "Sorry, file loading failed! Please check the file format.")    

### Table Events

def OnDoubleClick(event):
    item = table.selection()[0]
    value = table.item(item, 'values')
    geneID = value[0]
    uid = value[1]
    sequence = value[2]
    rc_sequence = reverseComplement(sequence)
    
    text_geneID.delete('1.0', tk.END)
    text_geneID.insert('1.0', str(geneID))
    
    text_uid.delete('1.0', tk.END)
    text_uid.insert('1.0', str(uid))
    
    text_sequence.delete('1.0', tk.END)
    text_sequence.insert('1.0', str(sequence))
    
    text_rc_sequence.delete('1.0', tk.END)
    text_rc_sequence.insert('1.0', str(rc_sequence))
    

def sortby(tree, col, descending):
    """sort tree contents when a column header is clicked on"""
    # grab values to sort
    data = [(tree.set(child, col), child) for child in tree.get_children('')]
    # if the data to be sorted is numeric change to float
    #data =  change_numeric(data)
    # now sort the data in place
    data.sort(reverse=descending)
    for ix, item in enumerate(data):
        tree.move(item[1], '', ix)
    # switch the heading so it will sort in the opposite direction
    tree.heading(col, command=lambda col=col: sortby(tree, col, int(not descending)))

def display_in_table():
    try:
        for a in df.index:
            row = df.ix[a]
            table.insert("", "end", "", values=tuple(row)) 
    except NameError:
        messagebox.showwarning("No Sequences to be Displayed", 
                               "Sorry, there's no loaded sequences to be displayed! Please load sequence file first.") 

### Other Button Functions

def clear():
    for i in table.get_children():
        table.delete(i)

def browse():
    start_time = time.time()
    clear()
    display_in_table()
    delta_time = time.time() - start_time
    
    text_time.delete('1.0', tk.END)
    text_time.insert('1.0', str(delta_time))           

def buttonExport():   
    if filenameSequences == '' or (filenameFASTQ == '' and len(filenameFASTQs) == 0):
        messagebox.showwarning("No File Loaded", 
                               "Sorry, no file loaded! Please choose sequence file and FASTQ file first.")
    else:
        try:
            len(df)
            len(reads)
            directory = filedialog.askdirectory()
            df.to_csv(directory + '/Counts of ' + filenameSequences.split('/')[-1] + ' matched with ' + filenameFASTQ.split('/')[-1] + '.csv')
            messagebox.showinfo("File Exported", "File of counted sequences successfully exported!")        
        except NameError:
            messagebox.showwarning("Error: No Counted DataFrame Generated", 
                               "Sorry, no effective counted DataFrame generated! Please check the previous workflow.")

def buttonAbout():
    about_root=tk.Tk()
    
    w = 380 # width for the Tk root
    h = 310 # height for the Tk root

    # get screen width and height
    ws = about_root.winfo_screenwidth() # width of the screen
    hs = about_root.winfo_screenheight() # height of the screen

    # calculate x and y coordinates for the Tk root window
    x = (ws/2) - (w/2)
    y = (hs/2) - (h/2)

    # set the dimensions of the screen 
    # and where it is placed
    about_root.geometry('%dx%d+%d+%d' % (w, h, x, y))
    about_root.title('About Bellfort Sequence Parser')  
    about_root.iconbitmap('dna.ico')

    label_author=tk.Label(about_root,text='Bellfort Sequence Parser Version 2.0', font=('tahoma', 9))
    label_author.place(x=90,y=30)

    label_author=tk.Label(about_root,text='Copyright (C) 2016', font=('tahoma', 9))
    label_author.place(x=125,y=60)
    
    label_author=tk.Label(about_root,text='Chen Lab', font=('tahoma', 9))
    label_author.place(x=150,y=90)
    
    label_author=tk.Label(about_root,text='Human Genome Sequencing Center', font=('tahoma', 9))
    label_author.place(x=80,y=120)
    
    label_author=tk.Label(about_root,text='Department of Molecular and Human Genetics', font=('tahoma', 9))
    label_author.place(x=50,y=150)
    
    label_author=tk.Label(about_root,text='Baylor College of Medicine', font=('tahoma', 9))
    label_author.place(x=110,y=180)
   

    button_okay=ttk.Button(about_root, width=15, text='OK', command=about_root.destroy)
    button_okay.place(x=130, y=235)

    about_root.mainloop()

### Batch Mode
#### FASTQ Files Loaded Button

def buttonBrowseFASTQs():
    global filenameFASTQs
    
    try:
        filenameFASTQs = filedialog.askopenfilenames(filetypes=(('FASTQ files', '*.fastq'), 
                                                              ('All files', '*.*')))
        text_fileFASTQ.delete('1.0', tk.END)
        text_fileFASTQ.insert('1.0', filenameFASTQ.split('/')[-1])
        
        # Reset the progress bar///////////////
        progressbar['value'] = 0
        progressbar_loadFASTQ['value'] = 0
        
        # Reset the percentage
        text_percentage.delete('1.0', tk.END)
        text_percentage.insert('1.0', str('0%'))
        
        indicator_preprocess = 0
    except:
        filenameFASTQs = ''

#### Batch Process Button Series

def loadFASTQ_batch(filenameFASTQ):
    global reads
    
    start_time = time.time() 
    
    f = open(filenameFASTQ)

    reads = []

    #try:
    while 1:
        name = f.readline().rstrip()
        sequence = f.readline().rstrip()
        f.readline()
        quality = f.readline().rstrip()

        if len(name) == 0:
            break

        union = name, sequence

        reads.append(union)

    f.close()

    end_time = time.time()
    delta_time = end_time - start_time

    text_time.delete('1.0', tk.END)
    text_time.insert('1.0', str(delta_time))  

    text_readNum.delete('1.0', tk.END)
    text_readNum.insert('1.0', str(len(reads)))  

    '''
    except:
        messagebox.showwarning("File Loading Failed", 
                               "Sorry, file loading failed! Please check the file format.")    '''

def preprocessFASTQ_batch():
    global reads, kmer_dict_reads, indicator_batch, gain_file
    
    try:  
        gotten = text_sequence_len.get('1.0', tk.END)
        k = int(gotten.rstrip())
        
        if k > len(reads[0][1]):
            messagebox.showwarning("Target Sequence Length Error", 
                                   "Sorry, the target sequence length is more than read length. Please check.")
        elif k < 3:
            messagebox.showwarning("Sequence Too Short", 
                                   "Sorry, the target sequence length is too short which will make the program running slowly. Please check.")
        elif filenameSequences == '':
            messagebox.showwarning("No Sequences Loaded", 
                                   "Sorry, no sequences loaded! Please load sequences first.")
        else:
            kmer_dict_reads = {}
                        
            gain = gain_file/(len(reads)*2)

            for read in reads:
                for i in range(len(read[1])-k+1):
                    kmer_dict_reads[read[1][i:i+k]] = set()
                    
                indicator_batch += gain 
                
            for read in reads:
                for i in range(len(read[1])-k+1):
                    kmer_dict_reads[read[1][i:i+k]].add(read)
                    
                indicator_batch += gain 

    except NameError:
        messagebox.showwarning("No FASTQ File Loaded", 
                               "Sorry, no loaded FASTQ file found! Please load FASTQ file first.")

def matchAll_batch():
    global  kmer_dict_reads, df, indicator_batch, gain
    
    try:       
        arr = np.array(df)

        for i in range(len(arr)):
            key1 = arr[i,2]
            key2 = reverseComplement(key1)

            try:
                n1 = len(kmer_dict_reads[key1])
            except KeyError:
                n1 = 0

            try:
                n2 = len(kmer_dict_reads[key2])
            except KeyError:
                n2 = 0

            arr[i, 4] += n1 + n2
            arr[i, 5] += 1


        df = pd.DataFrame(arr, columns = ['gene_id', 'UID', 'seq', 'Reserved', 'Count', 'Tag'])
        #df = df.set_index('UID', drop=False) 

    except NameError:
        messagebox.showwarning("No FASTQ Preprocessed or No Sequences Loaded", 
                               "Sorry, no FASTQ preprocess implemented or no sequences file loaded! Please preprocess FASTQ or load sequences first.")    

def batchProcess():
    global filenameFASTQs, indicator_batch, gain_file
    
    start_time = time.time() 
    
    if filenameFASTQs == '':
        messagebox.showwarning('No FASTQ Chosen', 
                               'Sorry, no FASTQ file chosen! Please browse and choose FASTQ file first.')   
        
    elif len(df) == 0:
        messagebox.showwarning("No Sequences Loaded", 
                                   "Sorry, no sequences loaded! Please load sequences first.")
    else:
        indicator_batch = 0
        gain_file = 100/len(filenameFASTQs)
        
        for filenameFASTQ in filenameFASTQs:
            loadFASTQ_batch(filenameFASTQ)
            preprocessFASTQ_batch()
            matchAll_batch()
            
        indicator_batch = 100
            
        messagebox.showinfo('Matching Completed', 
                                'Tada! Counting of sequences matched successfully completed!')  
        
    delta_time = time.time() - start_time
    
    text_time.delete('1.0', tk.END)
    text_time.insert('1.0', str(delta_time))   

def start_batch_thread(event):
    global batch_thread, indicator_batch
    batch_thread = threading.Thread(target=batchProcess)
    batch_thread.daemon = True
    
    progressbar_batch['value'] = indicator_batch
    text_percentage_batch.delete('1.0', tk.END)
    text_percentage_batch.insert('1.0', str(int(indicator_batch))+'%')
    
    batch_thread.start()
    root.after(20, check_batch_thread)

def check_batch_thread():
    if batch_thread.is_alive():
        progressbar_batch['value'] = indicator_batch
        text_percentage_batch.delete('1.0', tk.END)
        text_percentage_batch.insert('1.0', str(int(indicator_batch))+'%')
        
        root.after(20, check_batch_thread)    

## Main Flow

headers = ['gene_id', 'UID', 'seq', 'Reserved', 'count', 'tag']
header_widths = [280, 150, 350, 100, 80, 100]

root = tk.Tk()

indicator_preprocess = 0
indicator_loadSequences = 0
indicator_matchAll = 0
indicator_batch = 0
filenameSequences = ''
filenameFASTQ = ''
recordNum = 0
count = 0
df = pd.DataFrame([])

root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
#root.attributes('-fullscreen', True)
root.title('Bellfort Sequence Parser')
root.iconbitmap('dna.ico')


# Multicolumn Listbox/////////////////////////////////////////////////////////////////////////////
table = ttk.Treeview(height="20", columns=headers, selectmode="extended")
table.pack(padx=10, pady=20, ipadx=1200, ipady=100)

i = 1
for header in headers:
    table.heading('#'+str(i), text=header.title(), anchor=tk.W, command=lambda c=header: sortby(table, c, 0))
    table.column('#'+str(i), stretch=tk.NO, minwidth=0, width=tkf.Font().measure(header.title())+header_widths[i-1]) 
    i+=1    
table.column('#0', stretch=tk.NO, minwidth=0, width=0)

table.bind("<Double-1>", OnDoubleClick)
#///////////////////////////////////////////////////////////////////////////////////////////

# Scrollbar////////////////////////////////////////////////////////////////////////////////////////
vsb = ttk.Scrollbar(table, orient="vertical",  command = table.yview)
hsb = ttk.Scrollbar(table, orient="horizontal", command = table.xview)
## Link scrollbars activation to top-level object
table.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
## Link scrollbar also to every columns
map(lambda col: col.configure(yscrollcommand=vsb.set,xscrollcommand=hsb.set), table)
vsb.pack(side = tk.RIGHT, fill = tk.Y)
hsb.pack(side = tk.BOTTOM, fill = tk.X)        

#//////////////////////////////////////////////////////////////////////////////////////////////
y0 =310
y1 = 350
y2 = 420
y3 = 460
y4 = 520
y5 = 555
y6 = 595
y7 = 645
y8 = 695
# Text /////////////////////////////////////////////////////////////////////////////////////
text_recordNum=tk.Text(root, width=18, height=1, font=('tahoma', 9), bd=2, wrap='none')
text_recordNum.place(x=840, y=y0)
label_recordNum=tk.Label(root, text='records', font=('tahoma', 9))
label_recordNum.place(x=1000,y=y0)

text_fileSequences=tk.Text(root, width=55, height=1, font=('tahoma', 9), bd=2, wrap='none')
text_fileSequences.place(x=60, y=y0)

text_fileFASTQ=tk.Text(root, width=36, height=1, font=('tahoma', 9), bd=2, wrap='none')
text_fileFASTQ.place(x=60, y=y4)

text_count=tk.Text(root, width=16, height=1, font=('tahoma', 9), bd=2)
text_count.place(x=1000, y=y3)
label_count=tk.Label(root, text='Count:', font=('tahoma', 9))
label_count.place(x=940,y=y3)

text_geneID=tk.Text(root, width=20, height=1, font=('tahoma', 9), bd=2)
text_geneID.place(x=140, y=y2)
label_geneID=tk.Label(root, text='Gene ID:', font=('tahoma', 9))
label_geneID.place(x=60,y=y2)

text_uid=tk.Text(root, width=20, height=1, font=('tahoma', 9), bd=2)
text_uid.place(x=390, y=y2)
label_uid=tk.Label(root, text='UID:', font=('tahoma', 9))
label_uid.place(x=340,y=y2)

text_sequence=tk.Text(root, width=38, height=1, font=('tahoma', 9), bd=2)
text_sequence.place(x=680, y=y2)
label_sequence=tk.Label(root, text='Sequence:', font=('tahoma', 9))
label_sequence.place(x=600,y=y2)

text_rc_sequence=tk.Text(root, width=38, height=1, font=('tahoma', 9), bd=2)
text_rc_sequence.place(x=1000, y=y2)

text_sequence_len=tk.Text(root, width=5, height=1, font=('tahoma', 9), bd=2)
text_sequence_len.place(x=1260, y=y5)
label_sequence_len=tk.Label(root, text='nts', font=('tahoma', 9))
label_sequence_len.place(x=1315,y=y5)
text_sequence_len.delete('1.0', tk.END)
text_sequence_len.insert('1.0', str(20))

text_readNumDivided=tk.Text(root, width=13, height=1, font=('tahoma', 9), bd=2, wrap='none')
text_readNumDivided.place(x=335, y=y3+10)
label_readNumDivided1=tk.Label(root, text='by every', font=('tahoma', 9))
label_readNumDivided1.place(x=255,y=y3+10)
label_readNumDivided2=tk.Label(root, text='reads', font=('tahoma', 9))
label_readNumDivided2.place(x=460,y=y3+10)
text_readNumDivided.delete('1.0', tk.END)
text_readNumDivided.insert('1.0', str(250000))

text_readNum=tk.Text(root, width=22, height=1, font=('tahoma', 9), bd=2, wrap='none')
text_readNum.place(x=400, y=y6)
label_readNum=tk.Label(root, text='reads', font=('tahoma', 9))
label_readNum.place(x=590,y=y6)

text_time=tk.Text(root, width=15, height=1, font=('tahoma', 9), bd=2)
text_time.place(x=115, y=y8)
label_time=tk.Label(root, text='Time:', font=('tahoma', 9))
label_time.place(x=60,y=y8)
label_seconds=tk.Label(root, text='second(s)', font=('tahoma', 9))
label_seconds.place(x=250,y=y8)

text_percentage=tk.Text(root, width=8, height=1, font=('tahoma', 9), bg='gray95', bd=0)
text_percentage.place(x=1260, y=y4)

text_percentage_batch=tk.Text(root, width=8, height=1, font=('tahoma', 9), bg='gray95', bd=0)
text_percentage_batch.place(x=1260, y=y7)

# ProgressBar /////////////////////////////////////////////////////////////////////////////
progressbar_loadSequences = ttk.Progressbar(root, length=200, maximum=100, mode='determinate')
progressbar_loadSequences.place(x=530,y=y0)

progressbar_loadFASTQ = ttk.Progressbar(root, length=300, mode='indeterminate')
progressbar_loadFASTQ.place(x=400,y=y4)

progressbar = ttk.Progressbar(root, length=460, maximum=100, mode='determinate')
progressbar.place(x=760,y=y4)

progressbar_batch = ttk.Progressbar(root, length=520, maximum=100, mode='determinate')
progressbar_batch.place(x=700,y=y7)

# Button /////////////////////////////////////////////////////////////////////////////////
button_loadSequences = ttk.Button(root, text="Load sgRNA", width=20, command=loadSequences)
button_loadSequences.place(x=60, y=y1)

button_loadHalfMatchedSequences = ttk.Button(root, text="Load Half Matched sgRNA", width=30, command=buttonLoadHalfMatchedSequences)
button_loadHalfMatchedSequences.place(x=265, y=y1)

button_clear = ttk.Button(root, text="Clear", width=20, command=clear)
button_clear.place(x=1180, y=y1)

button_refresh = ttk.Button(root, text="Browse", width=20, command=browse)
button_refresh.place(x=1180, y=y0)

button_browseFASTQ = ttk.Button(root, text="Browse FASTQ...", width=20, command=buttonBrowseFASTQ)
button_browseFASTQ.place(x=60, y=y5)

button_divideFASTQ = ttk.Button(root, text="Divide FASTQ", width=20, command=lambda:start_divideFASTQ_thread(None))
button_divideFASTQ.place(x=60, y=y3+10)

button_loadFASTQ = ttk.Button(root, text="Load FASTQ", width=20, command=lambda:start_loadFASTQ_thread(None))
button_loadFASTQ.place(x=400, y=y5)

button_preprocessFASTQ = ttk.Button(root, text="Preprocess FASTQ & Count All Matched Sequences", width=55, command=lambda:start_preprocess_thread(None))
button_preprocessFASTQ.place(x=760, y=y5)

button_match = ttk.Button(root, text="Match", width=20, command=buttonMatch)
button_match.place(x=680, y=y3)

button_matchAll = ttk.Button(root, text="Match All", width=20, command=lambda:start_matchAll_thread(None))
button_matchAll.place(x=1180, y=y3)

button_about = ttk.Button(root, text="About", width=20, command=buttonAbout)
button_about.place(x=980, y=y8)

button_export = ttk.Button(root, text="Export", width=20, command=buttonExport)
button_export.place(x=720, y=y8)

button_exit = ttk.Button(root, text="Exit", width=20, command=root.destroy)
button_exit.place(x=1180, y=y8)

button_browseFASTQs = ttk.Button(root, text="Browse FASTQs...", width=25, command=buttonBrowseFASTQs)
button_browseFASTQs.place(x=60, y=y7)

button_batchProcess = ttk.Button(root, text="Batch Process", width=30, command=lambda:start_batch_thread(None))
button_batchProcess.place(x=400, y=y7)

root.bind('<Return>', start_preprocess_thread)
root.bind('<Return>', start_loadFASTQ_thread)
root.bind('<Return>', start_divideFASTQ_thread)
root.bind('<Return>', start_matchAll_thread)

root.mainloop()