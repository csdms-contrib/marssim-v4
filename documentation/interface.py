import sys
import os
#import shutil
from pathlib import Path
global SP
SP = os.sep
Firsttext = True
gdshow = ""
bibshow = ""
mssmshow = ""
nnn=1
mmm=1
bbb=1
global loc, pfiles, preffiles, docloc, bakloc
loc = "." + SP
pfiles = "parameter_files"
preffiles = "doc_file_prefixes"
docloc = loc + "doc"
bakloc = loc + "backup"
from appJar.appjar import gui
global app
def changeTab():
    pass
def do_documentation():
    pass
def do_afile(rb):
    """
    This reads and parses the lines in a parameter file
    """
    global nnn,mmm
    winname="YY"+str(nnn)
    panename="ZZ"+str(mmm)
    mmm=mmm+1
    nnn=nnn+1
    global thelist , xitems
    editfile = app.getRadioButton("Selectfile")
    xstart = 0
    for xstart in range(0,len(thefiles)):
        if editfile == thefiles[xstart]:
            theindex = xstart
    xstr = theprefixes[theindex]
    parfile = loc + SP + editfile
    firstup = True
    pindex = open(parfile,"r")
    # this writes the backup for the parameter file
    pbackup = loc + bakloc + SP + editfile + '.bak'
    pbackup_file = open(pbackup,"w")
    for param in pindex.readlines():
        pbackup_file.writelines(param)
    pbackup_file.close()
    pindex.close()
    pindex = open(parfile,"r")
    thelist = []
    maxvaluelen = 0
    for param in pindex.readlines():
        """
     This decodes the lines in the parameter file into the
     parameter values and the parameter explanations.
        """
        ptrim = param[0:-1]
        x = len(ptrim)
        comment_loc = ptrim.find("!")
        if comment_loc>0:
            the_value = ptrim[0:comment_loc]
            if len(the_value) > maxvaluelen:
                maxvaluelen=len(the_value)
            the_comment = ptrim[comment_loc+1:]
            if firstup:
                thelist.append ([the_value,the_comment])
                firstup = False
            else:
                thelist.append([the_value,the_comment])
        elif x>0:
            if firstup:
                thelist.append (['',ptrim])
                firstup = False
            else:
                thelist.append (['',ptrim])
                firstup = Falsethelist = []
    xitems = len(thelist)
    def exits():
        app.destroyScrollPane()
        app.destroySubWindow(winname)
    app.setBg('bisque')
    app.size = (1000,800)
    with app.subWindow("winname",title=editfile,blocking=True,transient=True, modal=True,grouped=True):
        """
        This is the routine to display the parameter editing file
        window
        """
        app.setSticky("nsew")
        app.size = (1000,600)
        app.setSticky("nsew")
        app.setBg('bisque')
        def dogo():
            return True
        def nogo():
            return False      
        app.startScrollPane(panename)
        app.size = (1000,600) 
        def doquit():
             app.destroySubWindow(thebut)
        def pressaction(but):
            """
            This reads and displays the explanatory text associated
            with pressing a parameter explanation button
            """
            thefile=loc + docloc + SP + but + '.txt'
            global thebut
            try:
                thetext =""
                handle = open(thefile,"r")
                thetext = handle.read()
                thebut=but +'x'
                thebut1= but + 'z'      
                with app.subWindow(thebut,blocking=True,transient=True, modal=True,grouped=True):
                    app.setSize(800,600)
                    app.setSticky("nsew")
                    app.startScrollPane("pressit")
                    app.setSubWindowBg=('gold')        
                    app.addMessage(thebut,thetext)
                    app.setMessageAspect(thebut,800)
                    app.setMessageBg=('gold')        
                    app.showSubWindow(thebut,hide=False)
                    app.destroySubWindow(thebut)
                    app.stopScrollPane()       
            except IOError:
                pass
        def dosave():
            """
            This saves the edited parameter file if the
            SAVE button is pressed
            """
            poutfile = open(parfile,"w")
            xcount = -1
            for x in thelist:
                xcount = xcount +1
                name = str(xcount)
                entr = 'entr'+name+'0'      
                if len(x[0])>0:        
                    parvalue = app.getEntry(entr)
                    parvalue= parvalue.ljust(maxvaluelen)
                else:
                    parvalue = ''
                if len(x[0] )>0:
                       theline = parvalue + "!"+ x[1] + "\n"
                else:
                       theline = x[1] + "\n"
                poutfile.writelines(theline)
            poutfile.close()
        xcount = -1
        for x in thelist:
            """ This displays the individual lines in the parameter
                File as the three columns: parameters, display buttion,
                and summary explanation
            """
            xcount=xcount+1
            for y in range(3):
                name = str(xcount)
                entr = 'entr'+name+str(y)
                but = xstr[0:2] + name+ 'x'
                lab = 'lab' + name +str(y)
                #print(entr)
                app.setSticky("w")
                if y == 0:
                    if len(x[0])>0:
                        app.addEntry(entr, row=xcount, column=1)
                        app.setEntry(entr,x[0],callFunction=False)
                if y == 1:
                    app.addButton(but,row=xcount,column=2,func=pressaction)
                    app.setButtonBg(but,'SkyBlue1')
                if y == 2:
                    app.addLabel(lab,text=x[1],row=xcount,column=3)
        app.setSticky("ew")
        app.addButton('Save', row=xcount+2,colspan=3,func=dosave)
        app.setButtonBg('Save','spring green')
        app.showSubWindow("winname")
        app.stopScrollPane()
        app.destroyAllSubWindows()
def do_edit():
    """
    This constructs the window containing the list of parameter files
    and displays the one that is selected for editing/display
    """
    global  rfile, xfile
    global thefiles , findex
    global theprefixes , xindex       
    app.addLabel("filex","SELECT THE PARAMETER FILE TO EDIT",row=3)
    rfile = docloc + SP + pfiles
    xfile = docloc + SP + preffiles
    #print (rfile)
    #print (xfile)
    findex = open(rfile,"r")
    xindex = open(xfile,"r")
    thefiles = []
    theprefixes = []
    while True:
        line = findex.readline()
        prefix =xindex.readline()
        if ("" == line):
            break
        else:
            stripline=line[0:-1]
            thefiles.append(stripline)
            theprefixes.append(prefix[0:-1])
    app.size = (800,600)
    app.stretch = 'both'
    app.setBg('bisque2')
    for x in range(0,len(thefiles)):
        app.addRadioButton("Selectfile",thefiles[x])
    app.setRadioButtonChangeFunction("Selectfile", do_afile)
def do_run():
    pass
def do_exit():
    pass
def head_val():
    """
    At the top of the parameter editing tab this will
    display the text discussing the editing procedure
    """
    global headitshow
    thedirectory=loc+docloc
    headittext=thedirectory + SP + "parameter_editing.txt"
    headfile=open(headittext,"r")
    headitshow=headfile.read() 
    with app.subWindow("headit",blocking=True,transient=True, modal=True,grouped=True,bg="plum1"):
        app.setFont(12)
        app.setBg('LightBlue1')        
        app.addMessage("headit",headitshow)       
        app.setMessageAspect("headit",350)
        app.showSubWindow("headit")
def write_entries():
    """
    This writes the default directories for retreiving and storing
    The parameter files
    """
    global direct,theinfile,theoutfile
    #print('writing with ',theinfile,theoutfile)
    thedirectory=loc+docloc
    direct= thedirectory + SP +"inout.txt"
    inoutfile = open(Path(direct),"w")
    toin=theinfile+'\n'
    toout=theoutfile+'\n'
    inoutfile.writelines(toin)
    inoutfile.writelines(toout)
    inoutfile.close()
def read_default_entries():
    """
    This reads the default directories for retreiving and storing
    The parameter files
    """
    global direct,theinfile,theoutfile
    thedirectory=loc+docloc
    direct= thedirectory + SP + "inout.txt"
    inoutfile = open(Path(direct),"r")
    theinfile=inoutfile.readline()
    theinfile=theinfile.replace("\n","")
    theoutfile=inoutfile.readline()
    theoutfile=theoutfile.replace("\n","")
    inoutfile.close()
    #print('reading with ',theinfile,theoutfile)    
def in_direct():
    """
    This transfers the parameter files from an existing directory to
    the editing directory.
    """
    global filedirectory,theinfile,theoutfile                 
    with app.subWindow("infile",blocking=True,transient=True,modal=True,grouped=True,bg="light pink"):
        def incopy_parameters():
            global filedirectory,theinfile,theoutfile
            theinfile=app.getEntry("iv")
            theinfile=theinfile.rstrip()
            isdir=os.path.isdir(Path(theinfile))
            if (isdir==True):
                thedirectory=loc
                rfile = docloc + SP + pfiles
                findex = open(rfile,"r")
                thefiles = []
                while True:
                    line = findex.readline()
                    if ("" == line):
                        break
                    else:
                        stripline=line[0:-1]
                        thefiles.append(stripline) 
                for x in range(0,len(thefiles)):
                   infile=theinfile + '/' + thefiles[x]
                   outfile=loc + SP + thefiles[x]
                   fromcopy=open(Path(infile),"r")
                   tocopy=open(outfile,"w")
                   for param in fromcopy.readlines():
                        tocopy.writelines(param)
                   fromcopy.close()
                   tocopy.close()
                write_entries()
            else:
                app.errorBox("FILE DIRECTORY ERROR","The input directory does not exist",parent="infile")
        app.setFont(12)
        read_default_entries()
        app.size=(600,200)
        app.setSticky("ew")
        app.addLabel("ix","Edit entry with / separator, no final separator")
        app.addEntry("iv")
        app.setEntry("iv",theinfile.ljust(80),callFunction=False)      
        app.addButton("Read Parameter Files",incopy_parameters)
        app.showSubWindow("infile")
        write_entries()
        app.destroySubWindow("infile")
def out_direct():
    """
    This transfers the parameter files from the editing directory to
    another directory, generally the execution director
    """
    global filedirectory,theinfile,theoutfile
    with app.subWindow("outfile",blocking=True,transient=True,modal=True,grouped=True,bg="light pink"):
        global filedirectory,theinfile,theoutfile
        def outcopy_parameters():
            global filedirectory,theinfile,theoutfile
            theoutfile=app.getEntry("ov")
            theoutfile=theoutfile.rstrip()
            isdir=os.path.isdir(Path(theoutfile))
            if (isdir==True):    
                thedirectory=loc
                rfile = docloc + SP + pfiles
                findex = open(rfile,"r")
                thefiles = []
                while True:
                    line = findex.readline()
                    if ("" == line):
                        break
                    else:
                        stripline=line[0:-1]
                        thefiles.append(stripline)
                for x in range(0,len(thefiles)):
                   infile=loc + SP + thefiles[x]
                   outfile=theoutfile +'/' + thefiles[x]
                   fromcopy=open(Path(infile),"r")
                   tocopy=open(outfile,"w")
                   for param in fromcopy.readlines():
                        tocopy.writelines(param)
                   fromcopy.close()
                   tocopy.close()
                write_entries()
            else:
                app.errorBox("FILE DIRECTORY ERROR","The output directory does not exist",parent="outfile")       
        app.setFont(12)
        read_default_entries()
        app.size=(600,100)
        app.setSticky("ew")        
        app.addLabel("ix","Edit Entry with / separator, no final separator")
        app.addEntry('ov')
        app.setEntry("ov",theoutfile.ljust(80),callFunction=False)        
        app.addButton("Write Parameter Files",outcopy_parameters)
        app.showSubWindow("outfile")        
        app.destroySubWindow("outfile")
def descr(which):
    """
    This controls display of the individual descriptive text
    windows of the leftmowt tab of the main interface window
    """
    global Firsttext, gdshow, bibshow,mssmshow
    if Firsttext:
        thedirectory=loc+docloc
        gdtext=thedirectory + SP + "gd.txt"
        bibtext=thedirectory + SP + "bib.txt"
        iotext =thedirectory + SP + "io.txt"
        mssmtext = thedirectory + SP + "mssm.txt"
        grdstext = thedirectory + SP + "grds.txt"
        statetxt = thedirectory + SP + "state_variables.txt"
        gdfile= open(gdtext,"r")
        gdshow = gdfile.read()
        bibfile = open(bibtext,"r")
        bibshow= bibfile.read()
        iofile= open(iotext,"r")
        ioshow = iofile.read()        
        mssmfile= open(mssmtext,"r")
        mssmshow = mssmfile.read()
        grdsfile = open(grdstext,"r")
        grdshow = grdsfile.read()
        statefile = open(statetxt,'r')
        stateshow = statefile.read()
    if which == "General Description":
        with app.subWindow("gd",blocking=True,transient=True, modal=True,grouped=True,bg='light goldenrod'):
            app.setFont(12)
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("gdpane")           
            app.addMessage("gd",gdshow)
            app.setMessageAspect("gd",400)
            app.showSubWindow("gd")
            app.stopScrollPane()            
    elif which == "Bibliography":
         with app.subWindow("bib",blocking=True,transient=True, modal=True,grouped=True,bg="powder blue"):
            app.setFont(12)
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("bibpane")            
            app.addMessage("bib",bibshow)
            app.setMessageAspect("bib",400)
            app.showSubWindow("bib")
            app.stopScrollPane()            
    elif which == "Input/Output Files":
        with app.subWindow("io",blocking=True,transient=True, modal=True,grouped=True,bg="pale green"):
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("iopane")
            app.setFont(12)            
            app.addLabel("io",ioshow)      
            app.showSubWindow("io")
            app.stopScrollPane()
    elif which == "Running MARSSIM":
         with app.subWindow("mssm",blocking=True,transient=True, modal=True,grouped=True,bg="plum1"):
            app.setFont(12)            
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("masmpane")
            app.addMessage("mssm",mssmshow)
            app.setMessageAspect("mssm",400)
            app.showSubWindow("mssm")
            app.stopScrollPane()            
    elif which == "Grids and Methods":
         with app.subWindow("grds",blocking=True,transient=True, modal=True,grouped=True,bg="plum1"):
            app.setFont(12)            
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("grdspane")
            app.addMessage("grds",grdshow)
            app.setMessageAspect("grds",400)
            app.showSubWindow("grds")
            app.stopScrollPane()           
    elif which == "State Variables":
         with app.subWindow("state",blocking=True,transient=True, modal=True,grouped=True,bg="plum1"):
            app.setFont(12)
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("statepane")
            app.addMessage("state",stateshow)
            app.setMessageAspect("state",400)
            app.showSubWindow("state")
            app.stopScrollPane()
    app.destroyAllSubWindows()
def showimage():
    """
    This is a stub for the popup that
    will show the current state of the simulation as a
    shaded relief image.  It will be part of the rightmost
    tab in the main GUI window
    """
    imagesizefile=loc + SP + "IMAGESIZE.DAT"
    sizefile = open (imagesizefile,"r")
    aline=''
    MX = 0
    MY = 0
    aline = sizefile.readline()
    aline = aline.strip()
    spaceloc = aline.find(" ")
    MX = int(aline[0:spaceloc])
    MY = int(aline[spaceloc+1:])
    print(MX,MY)
    with app.subWindow("ximage",blocking=True,transient=True, modal=True,grouped=True):
        app.size = (MX,MY)
        app.setSticky("nsew")        
        theimage = loc + SP + "B_pstag4280.pgm"
        app.addCanvas("OO")
        app.addCanvasImage("OO",x=MX//2, y=MY//2,image=theimage)
        app.showSubWindow("ximage")  
with gui("MARSSIM GUI",bg="blue") as app:
    """
    This is the main routine of the GUI, setting up the three
    Tabs and the contents contained within the tabs
    """
    app.setFont(family="Arial",weight="bold")
    with app.tabbedFrame("Tabs"):
        app.setTabbedFrameTabExpand("Tabs",expand=True)
        app.size = (800,300)
        with app.tab("DOCUMENTATION",bg="sky blue"):
            app.addLabel("l1","MARSSIM is a landform evolution model that is oriented towards simulating multiple planetary processes.",colspan=2)
            app.button("General Description",descr,pos=(1,0))
            app.button("Bibliography",descr,pos=(1,1))
            app.button("Input/Output Files",descr,pos=(2,0))
            app.button("Running MARSSIM",descr,pos=(2,1))
            app.button("Grids and Methods",descr,pos=(3,0))
            app.button("State Variables",descr,pos=(3,1))
        with app.tab("EDIT PARAMETERS",bg="goldenrod"):
            app.setSticky("ew")
            app.button("Parameter Editing and Description",head_val,row=0,colspan=2,bg="coral")
            app.button("Read parameter files from a directory",in_direct,pos=(1,0),bg='salmon')
            app.button("Write parameter files to a directory",out_direct,pos=(1,1),bg='green2')        
            do_edit()
        #with app.tab("RUN MARSSIM",bg="sienna1"):
            #app.addButton("Show Current State",func=showimage)
            #do_run()

