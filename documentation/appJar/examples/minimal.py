import sys
import time
sys.path.append("../../")
Firsttext = True
gdshow = ""
bibshow = ""
mssmshow = ""
nnn=1
mmm=1
bbb=1
global loc, pfiles, preffiles, docloc, bakloc
loc = r"E:\python\marssim\parameters" + "\\"
pfiles = "parameter_files"
preffiles = "doc_file_prefixes"
docloc = r"\doc"
bakloc = r"\backup"
from appJar import gui
global app
def changeTab():
    pass
def do_documentation():
    pass
def do_afile(rb):
    global nnn,mmm
    winname="YY"+str(nnn)
    panename="ZZ"+str(mmm)
    mmm=mmm+1
    nnn=nnn+1
    global thelist , xitems
    editfile = app.getRadioButton("Selectfile")
    print(editfile)
    xstart = 0
    for xstart in range(0,len(thefiles)):
        if editfile == thefiles[xstart]:
            theindex = xstart
            #print  (theindex)
    xstr = theprefixes[theindex]
    parfile = loc + "\\" + editfile
    print (parfile)
    firstup = True
    pindex = open(parfile,"r")
    pbackup = loc + bakloc + "\\" + editfile + '.bak'
    pbackup_file = open(pbackup,"w")
    for param in pindex.readlines():
        pbackup_file.writelines(param)
    pbackup_file.close()
    pindex.close()
    pindex = open(parfile,"r")
    thelist = []
    maxvaluelen = 0
    for param in pindex.readlines():
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
        app.setSticky("nsew")
        app.size = (1000,600)
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
            thefile=loc + docloc + "\\" + but + '.txt'
            global thebut
            try:
                thetext =""
                handle = open(thefile,"r")
                thetext = handle.read()
                thebut=but +'x'
                thebut1= but + 'z'      
                with app.subWindow(thebut,blocking=True,transient=True, modal=True,grouped=True):
                    app.setSize(800,600)
                    app.setSubWindowBg=('gold')        
                    app.addMessage(thebut,thetext)
                    app.setMessageAspect(thebut,800)
                    app.setMessageBg=('gold')        
                    app.showSubWindow(thebut,hide=False)                    
                    app.stopSubWindow()        
            except IOError:
                pass
        def dosave():
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
            xcount=xcount+1
            for y in range(3):
                name = str(xcount)
                entr = 'entr'+name+str(y)
                but = xstr[0:2] + name+ 'x'
                lab = 'lab' + name +str(y)
                #print(entr)
                if y == 0:
                    if len(x[0])>0:
                        app.addEntry(entr, row=xcount, column=1)
                        app.setEntry(entr,x[0],callFunction=False)
                if y == 1:
                    app.addButton(but,row=xcount,column=2,func=pressaction)
                    app.setButtonBg(but,'SkyBlue1')
                if y == 2:
                    app.addLabel(lab,text=x[1],row=xcount,column=3)
        app.addButton('Save', row=xcount+2,column=1,func=dosave)
        app.setButtonBg('Save','spring green')
        app.showSubWindow("winname")
        app.stopScrollPane()
        app.destroyAllSubWindows()
def do_edit():
    global  rfile, xfile
    global thefiles , findex
    global theprefixes , xindex       
    xloc = app.getEntry("Directory")
    app.addLabel("filex","SELECT THE PARAMETER FILE TO EDIT",row=3)
    rfile = loc + docloc + "\\" + pfiles
    xfile = loc + docloc + "\\" + preffiles
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
    #print (thefiles)
    pass
def do_run():
    pass
def do_exit():
    pass
def descr(which):
    global Firsttext, gdshow, bibshow,mssmshow
    if Firsttext:
        thedirectory=loc+docloc
        gdtext=thedirectory + "\\" + "gd.txt"
        bibtext=thedirectory + "\\" + "bib.txt"
        iotext =thedirectory + "\\" + "io.txt"
        mssmtext = thedirectory + "\\" + "mssm.txt"
        grdstext = thedirectory + "\\" + "grds.txt"
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
    if which == "General Description":
        with app.subWindow("gd",blocking=True,transient=True, modal=True,grouped=True,bg='light goldenrod'):
            app.setFont(12)            
            app.addMessage("gd",gdshow)
            app.setMessageAspect("gd",400)
            app.showSubWindow("gd") 
    elif which == "Bibliography":
         with app.subWindow("bib",blocking=True,transient=True, modal=True,grouped=True,bg="powder blue"):
            app.setFont(12)            
            app.addMessage("bib",bibshow)
            app.setMessageAspect("bib",400)
            app.showSubWindow("bib")
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
            app.addMessage("mssm",mssmshow)
            app.setMessageAspect("mssm",400)
            app.showSubWindow("mssm") 
    elif which == "Grids and Methods":
         with app.subWindow("grds",blocking=True,transient=True, modal=True,grouped=True,bg="plum1"):
            app.setFont(12)            
            app.addMessage("grds",grdshow)
            app.setMessageAspect("grds",400)
            app.showSubWindow("grds")
    app.destroyAllSubWindows()
with gui("MARSSIM GUI",bg="blue") as app:
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
        #app.stopTab()
        with app.tab("EDIT PARAMETERS",bg="goldenrod"):
            app.addLabel("Select the Directory containing the parameter files",row=0)
            app.addEntry("Directory",row=1)
            app.setEntrySubmitFunction("Directory",do_edit)
        with app.tab("RUN MARSSIM",bg="sienna1"):
            app.addLabel("l3","RUN MARSSIM")
            do_run()
        with app.tab("EXIT",bg="medium orchid"):
            app.addLabel("l4","EXIT")
            do_exit()
