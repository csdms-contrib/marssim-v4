import sys
import time
sys.path.append("../../")
Firsttext = True
gdshow = ""
bibshow = ""
mssmshow = ""
from appJar import gui
def changeTab():
    pass
def do_documentation():
    pass
def do_afile(rb):
    editfile = app.getRadioButton("Selectfile")
    print(editfile)
    xstart = 0
    for xstart in range(0,len(thefiles)):
        if editfile == thefiles[xstart]:
            theindex = xstart
            print  (theindex)
    
def do_edit():
    loc = app.getEntry("Directory")
    app.addLabel("filex","SELECT THE PARAMETER FILE TO EDIT",row=3)
    loc = r"E:\python\marssim\parameters" + "\\"
    pfiles = "parameter_files"
    preffiles = "doc_file_prefixes"
    docloc = r"\doc"
    bakloc = r"\backup"
    rfile = loc + docloc + "\\" + pfiles
    xfile = loc + docloc + "\\" + preffiles
    #rfile = loc + pfiles
    #xfile = loc + preffiles
    findex = open(rfile,"r")
    xindex = open(xfile,"r")
    global thefiles
    global theprefixes
    thefiles = []
    theprefixes = []
    while True:
        line = findex.readline()
        prefix =xindex.readline()
        if ("" == line):
            break
        else:
            stripline=line[0:-1]
            #print(stripline)
            thefiles.append(stripline)
            theprefixes.append(prefix[0:-1])
    #print(thefiles)
    #print (len(thefiles))
    #print (theprefixes)
    #print (len(theprefixes))
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
        Firsttext = False
        gdfile= open(r"gd.txt","r")
        gdshow = gdfile.read()
        bibfile = open(r"bib.txt","r")
        bibshow= bibfile.read()
        iofile= open(r"io.txt","r")
        ioshow = iofile.read()        
        mssmfile= open(r"mssm.txt","r")
        mssmshow = mssmfile.read()               
    if which == "General Description":
        def gdquit():
            app.destroyAllSubWindows() 
        with app.subWindow("gd",blocking=True,transient=False, modal=True,grouped=True,bg='light goldenrod'):
            app.setFont(12)            
            app.addMessage("gd",gdshow)
            app.setMessageAspect("gd",400)
            
            app.button("QUIT",gdquit)
            app.showSubWindow("gd",hide=False)
            #app.hideSubWindow("gd")
            app.stopSubWindow()
            #app.destroySubWindow("gd")
            app.go()
        print("gd")      
    elif which == "Bibliography":
         with app.subWindow("bib",blocking=True,transient=False, modal=True,grouped=True,bg="powder blue"):
            app.setFont(12)            
            app.addMessage("bib",bibshow)
            app.setMessageAspect("bib",400)
            def gdquit():
                app.destroySubWindow("bib")
            app.button("QUIT",gdquit)
            app.showSubWindow("bib")
            #app.hideSubWindow("bib")
            app.stopSubWindow()
            #app.destroySubWindow("bib")
            app.go()       
        #print("bib")
    elif which == "Input/Output Files":
        with app.subWindow("io",blocking=True,transient=False, modal=True,grouped=True,bg="pale green"):
            app.size = (1000,600)
            app.setSticky("nsew")
            app.startScrollPane("iopane")
            app.ScrollPane(disabled="horizontal")
            #app.width = 800
            #app.SrollPaneSize = (800,600)
            app.setFont(12)            
            app.addLabel("io",ioshow)
 #           app.startScrollPane("iopane",column=800,sticky="ew")
 #           app.ScrollPane(disabled="horizontal")           
            #app.setMessageAspect("io",400)
            #def ioquit():
                #app.destroySubWindow("io")
            #app.button("QUIT",ioquit)
            app.showSubWindow("io")
            #app.hideSubWindow("io")
            #app.stopSubWindow()
            #app.destroySubWindow("io")
            #app.stopSubWindow("io")
            app.stopScrollPane()
            app.go()
        #print("I//O")
    elif which == "Running MARSSIM":
         with app.subWindow("mssm",blocking=True,transient=False, modal=True,grouped=True,bg="plum1"):
            app.setFont(12)            
            app.addMessage("mssm",mssmshow)
            app.setMessageAspect("mssm",400)
            def gdquit():
                app.destroySubWindow("mssm")
            app.button("QUIT",gdquit)
            app.showSubWindow("mssm")
            #app.hideSubWindow("bib")
            app.stopSubWindow()
            #app.destroySubWindow("bib")
            app.go()       
        #print("bib")        
        #print("run")
    if which == "Grids and Methods":   
        print("grd")
    app.destroyAllSubWindows()
with gui("MARSSIM GUI",bg="blue") as app:
    #app.size = (800,500)
    #app.addMenuList("Tabs", ["EXIT", "DOCUMENTATION", "EDIT PARAMA", "RUN MARSSIM"], changeTab)
    #app.setTabbedFrameDisableAllTabs("Tabs", False)
    with app.tabbedFrame("Tabs"):
        app.setTabbedFrameTabExpand("Tabs",expand=True)
        with app.tab("DOCUMENTATION",bg="sky blue"):
            #,beforeTab=None,afterTab="EDIT PARAMETERS")
            app.addLabel("l1","MARSSIM is a landform evolution model that is oriented towards simulating multiple planetary processes.",colspan=2)
            app.button("General Description",descr,pos=(1,0))
            app.button("Bibliography",descr,pos=(1,1))
            app.button("Input/Output Files",descr,pos=(2,0))
            app.button("Running MARSSIM",descr,pos=(2,1))
            app.button("Grids and Methods",descr,pos=(3,0))
            #          
            #do_documentation()
        #app.stopTab()
        with app.tab("EDIT PARAMETERS",bg="goldenrod"):
            #,beforeTab="DOCUMENTATION",afterTab="RUN MARSSIM")
            #app.addLabel("l2","EDIT PARAMETERS")
            app.addLabel("Select the Directory containing the parameter files",row=0)
            app.addEntry("Directory",row=1)
            app.setEntrySubmitFunction("Directory",do_edit)
            #app.setSubmitFunction("Directory",do_edit)
            #loc = app.getEntry("Directory")
            #pfiles = "parameter_files"
            #preffiles = "doc_file_prefixes"
            #rfile = loc + '\\' + pfiles
            #xfile = loc + '\\' + preffiles
            #findex = open(rfile,"r")
            #xindex = open(xfile,"r")
            #thefiles = findex.readlines()
            #theprefixes = xindex.readlines()
            #nnn = 0
            #for lines in findex.readlines():
                #pfile = lines[0:-1]
                #nnn = nnn +1
                #xpref(lines) = xindex[0:-1]
                #print (pfile)
            #do_edit()
        #app.stopTab()
        with app.tab("RUN MARSSIM",bg="sienna1"):
            #,beforeTab="EDIT PARAMETERS",afterTab="EXIT")
            app.addLabel("l3","RUN MARSSIM")
            do_run()
        #app.stopTab()
        with app.tab("EXIT",bg="medium orchid"):
            #,beforeTab="RUN MARSSIM",afterTab=None)
            app.addLabel("l4","EXIT")
            do_exit()
        #app.stopTab()
        #pass
        #app.stopTabbedFrame()        
    #app.addToolbar(["EXIT", "DOCUMENTATION", "EDIT PARAMS","RUN MARSSIM"], toolbar, findIcon=False)
    #app.addMenuWindow()
    #app.addMenuHelp(toolbar)
    pass
    #logout()

#app.go()
