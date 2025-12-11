import ROOT
import sys
from array import array
import pickle
import ctypes
import math
import common.tdrstyle as tdrstyle
import common.CMS_lumi as CMS_lumi
from python.VHTools.py_helpers import *
sty = tdrstyle.setTDRStyle()

### = commented out to test tdrStyle, uncomment if not using

drawprelim = True

@contextmanager
def open_root_file(filename, mode="READ"):
    """
    Context manager for safely opening and closing a ROOT TFile.

    This ensures that the ROOT file is automatically closed when the
    context exits, even in the event of an exception.

    Parameters:
    ----------
    filename : str
        Path to the ROOT file to open.

    mode : str, optional
        Mode in which to open the file (default is "READ").
        Typical options include:
        - "READ"   : Open existing file for reading
        - "UPDATE" : Open existing file for reading/writing
        - "RECREATE": Create new file, overwriting if exists
        - "NEW"    : Create new file, fail if exists

    Yields:
    ------
    ROOT.TFile
        The opened ROOT file object.

    Example:
    --------
    >>> with open_root_file("data.root") as f:
    ...     h = f.Get("myHistogram")
    ...     h.Draw()
    """
    file = ROOT.TFile.Open(filename, mode)
    try:
        yield file
    finally:
        if file and file.IsOpen():
            file.Close()

class plotter_base(object):

    def __init__(self):
        self.fillstyle=1001
        self.linestyle=1
        self.linecolor=1
        self.linewidth=2
        self.fillcolor=ROOT.kOrange-3
        self.markerstyle=20
        self.corrFactors=[]

    def addCorrectionFactor(self,value,model):
        corr=dict()
        corr['value']=value
        corr['model']=model
        self.corrFactors.append(corr)

    def setLineProperties(self,linestyle,linecolor,linewidth):
        self.linestyle=linestyle
        self.linecolor=linecolor
        self.linewidth=linewidth 

    def setFillProperties(self,fillstyle,fillcolor):
        self.fillstyle=fillstyle
        self.fillcolor=fillcolor

    def setMarkerProperties(self,markerstyle):
        self.markerstyle=markerstyle



class rdf_plotter(plotter_base):
    def __init__(self,file,isMC=False,weight = "1.0",tree='Events',defaultCuts = "1.0",report='report'):
        self.rdf = ROOT.RDataFrame(tree,file)
        self.weight=weight
        self.isMC=isMC
        self.defaultCuts = defaultCuts
        self.file = file
        self.report = report
        super(rdf_plotter,self).__init__()
        #if MC read the sum of weights and weigh the event
        if isMC:
            with open_root_file(file) as f:
                t = f.Get('Runs')
                sumw = 0.0
                for event in t:
                    sumw += event.genEventSumw
                self.weight = self.weight+'*(genWeight/{})*(sample_sigma)*(Pileup_weight)'.format(sumw if sumw>0 else 1.0)

    def readReport(self) -> dict[str]:
        out = {}
        with open_root_file(self.file) as f:
            report = f.Get(self.report)
            if not report:
                raise RuntimeError(f"Report {self.report} not found in {self.file}")
            
            cuts = [b.GetName() for b in report.GetListOfBranches()]
            out = {c: 0 for c in cuts}
            for event in report:
                for c in cuts:
                    out[c] += getattr(event, c, 0)

            return out
        

    def define(self, var, definition):
        self.rdf = self.rdf.Define(var, definition)

    def redefine(self, var, definition):
        self.rdf = self.rdf.Redefine(var, definition)
    
    def display(self,var):
        self.rdf.Display(var).Print()

    def filter(self,condition):
        '''
        Inputs
        ======
        condition(string): condition to be used for filtering the dataframe that evaluates to True (1) or False (0)
        new_var(string): name of the new column  
        
        Returns
        =======
        RDataFrame with filter applied
        '''
        self.rdf = self.rdf.Filter(condition)

        
    def hist1d(self,var,cuts,lumi,model,titlex = "",units = ""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+lumi+"*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo1D(model,var,'plot_weight')
        h.Sumw2()
        #h.Sumw2(0)
        #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        for n in range(1, h.GetNbinsX()+1):
            if h.GetBinContent(n) < 0:
                h.SetBinContent(n, 0)
        return h


    def profile1d(self,var1,var2,cuts,lumi,model,titlex = "",unitsx = "",titley = "",unitsy = ""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+lumi+"*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Profile1D(model,var1,var2,'plot_weight')
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if unitsx=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        if unitsy=="":
            h.GetYaxis().SetTitle(titley)
        else:
            h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h


    def hist2d(self,var1,var2,cuts,lumi,model,titlex = "",unitsx = "",titley="",unitsy=""):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+lumi+"*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo2D(model,var1,var2,'plot_weight')
        h.Sumw2()
        #h.Sumw2(0)
        #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if unitsx=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
        if unitsy=="":
            h.GetYaxis().SetTitle(titley)
        else:
            h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
        return h

    def unrolled2d(self,var1,var2,cuts,lumi,model):
        corrString="1.0"
        for corr in self.corrFactors:
            corrString = corrString+"*("+str(corr['value'])+")" 
        c = "("+self.defaultCuts+")*("+cuts+")*"+lumi+"*"+self.weight+"*("+corrString+")"
        rdf=self.rdf.Define('plot_weight',c)
        h=rdf.Histo2D(model,var1,var2,'plot_weight')
        h.Sumw2()        
        hu=ROOT.TH1D(h.GetName(),h.GetTitle(),h.GetNbinsX()*h.GetNbinsY(),0,h.GetNbinsX()*h.GetNbinsY())
        for i in range(1,h.GetNbinsX()+1):
            for j in range(1,h.GetNbinsY()+1):
                b=h.GetBin(i,j)
                bin1d=(i-1)*h.GetNbinsX()+j
                hu.SetBinContent(bin1d,h.GetBinContent(b))
                hu.SetBinError(bin1d,h.GetBinError(b))
            
        hu.SetLineStyle(self.linestyle)
        hu.SetLineColor(self.linecolor)
        hu.SetLineWidth(self.linewidth)
        hu.SetFillStyle(self.fillstyle)
        hu.SetFillColor(self.fillcolor)
        hu.SetMarkerStyle(self.markerstyle)
        return hu


    

class merged_plotter(plotter_base):

    def __init__(self, plotters):
        self.fillstyle=1001
        self.linestyle=1
        self.linecolor=1
        self.linewidth=2
        self.fillcolor=ROOT.kOrange-3
        self.markerstyle=20
        self.corrFactors=[]
        self.plotters = plotters

    def define(self, var, definition):
        for plotter in self.plotters:
            plotter.define(var, definition)

    def filter(self, condition):
        for plotter in self.plotters:
            plotter.filter(condition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter.redefine(var, definition)

    def hist1d(self,var,cuts,lumi,model,titlex = "",units = ""):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.hist1d(var, cuts, lumi, model, titlex, units)
            else:
                h.Add(plotter.hist1d(var, cuts, lumi, model, titlex, units).GetValue())
        if h is None:
            tmprdf = ROOT.RDataFrame(1)
            tmprdf = tmprdf.Define("x", "0")
            h = tmprdf.Define("goodX", "x!=0").Histo1D(model, "goodX")
        h.Sumw2()
        #h.Sumw2(0)
        #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
        h.SetLineStyle(self.linestyle)
        h.SetLineColor(self.linecolor)
        h.SetLineWidth(self.linewidth)
        h.SetFillStyle(self.fillstyle)
        h.SetFillColor(self.fillcolor)
        h.SetMarkerStyle(self.markerstyle)
        if units=="":
            h.GetXaxis().SetTitle(titlex)
        else:
            h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
        return h


    def hist2d(self,var1,var2,cuts,lumi,model,titlex = "",unitsx = "",titley="",unitsy=""):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.hist2d(var1,var2, cuts, lumi, model, titlex, unitsx, titley, unitsy)
            else:
                h.Add(plotter.hist2d(var1,var2, cuts, lumi, model, titlex, unitsx, titley, unitsy).GetValue())
        if h is None:
            return h
        else:
            h.Sumw2()
            #h.Sumw2(0)
            #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
            h.SetLineStyle(self.linestyle)
            h.SetLineColor(self.linecolor)
            h.SetLineWidth(self.linewidth)
            h.SetFillStyle(self.fillstyle)
            h.SetFillColor(self.fillcolor)
            h.SetMarkerStyle(self.markerstyle)
            if unitsx=="":
                h.GetXaxis().SetTitle(titlex)
            else:
                h.GetXaxis().SetTitle(titlex+ " ["+unitsx+"]")
            if unitsy=="":
                h.GetYaxis().SetTitle(titley)
            else:
                h.GetYaxis().SetTitle(titley+ " ["+unitsy+"]")
            return h

    def unrolled2d(self,var1,var2,cuts,lumi,model):
        h = None
        for plotter in self.plotters:
            if h is None:
                h = plotter.unrolled2d(var1,var2, cuts, lumi, model)
            else:
                h.Add(plotter.unrolled2d(var1,var2, cuts, lumi, model))
        if h is None:
            return h
        else:
            h.Sumw2()
            #h.Sumw2(0)
            #h.SetBinErrorOption(ROOT.TH1D.kPoisson)
            h.SetLineStyle(self.linestyle)
            h.SetLineColor(self.linecolor)
            h.SetLineWidth(self.linewidth)
            h.SetFillStyle(self.fillstyle)
            h.SetFillColor(self.fillcolor)
            h.SetMarkerStyle(self.markerstyle)
        return h
                               
    def getArray(self, var, cuts, lumi):
        """
        Applies cuts to RDF and extracts numpy array of values for a specific variable.

        Parameters:
            vars (str): variable
            cuts (str): Selection cuts
            lumi (str): Luminosity weight
            range (tuple) : range over

        Returns:
            list of values corresponding to cuts
        """        
        # Initialize empty lists to store the data
        combined_var = []

        # Iterate over the plotters and collect data
        for plotter in self.plotters:
            # Define weight expression
            corrString = "1.0"
            for corr in plotter.corrFactors:
                corrString += f"*({corr['value']})"

            # Apply cuts and filters based on the provided cuts and lumi
            weight_expr = f"({plotter.defaultCuts})*({cuts})*{lumi}*{plotter.weight}*({corrString})"
            rdf = plotter.rdf.Define("plot_weight", weight_expr)
            rdf = rdf.Filter("plot_weight")

            # Get data from the current plotter
            numpy_var = rdf.AsNumpy([var])[var]

            # Collect the data
            combined_var.extend(numpy_var)

        # If no data points were collected, return None
        if len(combined_var) == 0:
            return None

        

        return combined_var
    
    def display(self,var):
        for plotter in self.plotters:
            plotter.display(var)

    def readReport(self):
        out = {}
        for plotter in self.plotters:
            r = plotter.readReport()
            for cut in r:
                if cut in out:
                    out[cut] += r[cut]
                else:
                    out[cut] = r[cut]
        return out



#Tricky, use with care!
class background_plotter(merged_plotter):
    def __init__(self,cutsSR,cutsCR,cutsSB,plotters):
        self.cutsSR=cutsSR
        self.cutsCR=cutsCR
        self.cutsSB=cutsSB
        super(background_plotter,self).__init__(plotters)
        self.define('scaleVar','1.0')
    def getScale(self,cuts,lumi='1.0'):       
        hNum=super(background_plotter,self).hist1d('scaleVar','&&'.join([cuts,self.cutsSB]),lumi,('a','a',10,-1,1))
        #replace cuts with cuts from control region
        hDenom=super(background_plotter,self).hist1d('scaleVar','&&'.join([cuts.replace(self.cutsSR,self.cutsCR),self.cutsSB]),lumi,('a','a',10,-1,1))

        if hDenom.Integral()==0:
            print("Error , denominator has zero events, returning scale=1")
            return 1.0
        else:
            return hNum.Integral()/hDenom.Integral()        
        
    def hist1d(self,var,cuts,lumi,model,titlex = "",units = ""):
        h=super(background_plotter,self).hist1d(var,cuts.replace(self.cutsSR,self.cutsCR),lumi,model,titlex,units)
        h.Scale(self.getScale(cuts,lumi))
        return h
    
    def hist2d(self,var1,var2,cuts,lumi,model,titlex = "",unitsx = "",titley="",unitsy=""):
        h=super(background_plotter,self).hist2d(var1,var2,cuts.replace(self.cutsSR,self.cutsCR),lumi,model,titlex,unitsx,titley,unitsy)
        h.Scale(self.getScale(cuts,lumi))
        return h

    def unrolled2d(self,var1,var2,cuts,lumi,model):
        h=super(background_plotter,self).unrolled2d(var1,var2,cuts.replace(self.cutsSR,self.cutsCR),lumi,model)      
        h.Scale(self.getScale(cuts,lumi))
        return h



    
class combined_plotter(object):
    def __init__(self,defaultCut="1"):
        self.plotters = []
        self.types    = []
        self.labels   = []
        self.names    = []
        self.log=False
        self.defaultCut=defaultCut

    def setLog(self,doLog):
        self.log=doLog

    def add_plotter(self,plotter,name="",label = "label",typeP = "background"):
        self.plotters.append(plotter)
        self.types.append(typeP)
        self.labels.append(label)
        self.names.append(name)

    def define(self, var, definition):
        for plotter in self.plotters:
            plotter.define(var, definition)

    def redefine(self, var, definition):
        for plotter in self.plotters:
            plotter.redefine(var, definition)


    def draw_stack(self,var,cut,lumi,model, titlex = "", units = "",expandY=0.0,SFs="(1)", verbose = False, prelim = "Work in progress", lumi_label = "", outOfFrame = 1):
###        canvas = ROOT.TCanvas("canvas","")
        canvas = ROOT.TCanvas("canvas", "", 800, 600)
#        ROOT.gStyle.SetOptStat(0)
#        ROOT.gStyle.SetOptTitle(0)
#        canvas.Range(-68.75,-7.5,856.25,42.5)
#        canvas.SetFillColor(0)
#        canvas.SetBorderMode(0)
#        canvas.SetBorderSize(2)
#        canvas.SetTickx(1)
#        canvas.SetTicky(1)
#        canvas.SetLeftMargin(0.15)
###        canvas.SetRightMargin(0.05) #commented to test tdrStyle
#        canvas.SetTopMargin(0.05)
#        canvas.SetBottomMargin(0.15)
#        canvas.SetFrameFillStyle(0)
#        canvas.SetFrameBorderMode(0)
#        canvas.SetFrameFillStyle(0)
#        canvas.SetFrameBorderMode(0)

        canvas.SetTopMargin(0.055) # testing tdr style
        canvas.SetBottomMargin(0.12) # testing
        canvas.SetLeftMargin(0.12) # testing

        canvas.cd()
        hists=[]
        stack = ROOT.THStack("stack","")
        
        signal=0
        background=0
        backgroundErr=0
        
        data=None
        dataG=None
        error=ctypes.c_double(0.0)

        cutL="("+self.defaultCut+")*("+cut+")"

        for (plotter,typeP,label,name) in zip(self.plotters,self.types,self.labels,self.names):
            if typeP == "signal" or typeP =="background":
                if plotter.isMC:
                    hist = plotter.hist1d(var,cutL+"*("+SFs+")",lumi,model,titlex,units)
                    hist.SetName(name)
                    stack.Add(hist.GetValue())
                    hists.append(hist.GetValue())
                else:
                    hist = plotter.hist1d(var,cutL+"*("+SFs+")",'1',model,titlex,units)
                    hist.SetName(name)
                    stack.Add(hist.GetValue())
                    hists.append(hist.GetValue())
                if verbose:
                    print( label+" : %f\n" % hist.Integral())
 
                if typeP == "signal" :
                    signal+=hist.Integral()
                if typeP == "background" :
                    background+=hist.IntegralAndError(1,hist.GetNbinsX(),error)
                    backgroundErr+=error.value*error.value
       
            if typeP =="data":
                hist = plotter.hist1d(var,cutL,"1",model,titlex,units)
                hist.SetName(hist.GetName()+label)
                hists.append(hist.GetValue())
                data=hist.GetValue()
                dataG=convertToPoisson(hist.GetValue())
                dataG.SetLineWidth(1)
                if verbose:
                    print( label+" : %f\n" % hist.Integral())
       
       
        #if data not found plot stack only
        if data != None:                  
            datamax = ROOT.Math.chisquared_quantile_c((1-0.6827)/2.,2*(data.GetMaximum()+1))/2.

        else: 
            datamax = stack.GetMaximum()
        if not self.log:
            frame = canvas.DrawFrame(hists[0].GetXaxis().GetXmin(),0.0,hists[0].GetXaxis().GetXmax(),max(stack.GetMaximum(),datamax)*(1.20+expandY*0.3))
        else:    
            frame = canvas.DrawFrame(hists[0].GetXaxis().GetXmin(),0.1,hists[0].GetXaxis().GetXmax(),max(stack.GetMaximum(),datamax)*100)

#        frame.GetXaxis().SetLabelFont(42)
#        frame.GetXaxis().SetLabelOffset(0.007)
#        frame.GetXaxis().SetLabelSize(0.045)
        frame.GetXaxis().SetTitleSize(0.05)
#        frame.GetXaxis().SetTitleOffset(1.15)
#        frame.GetXaxis().SetTitleFont(42)
#        frame.GetYaxis().SetLabelFont(42)
#        frame.GetYaxis().SetLabelOffset(0.007)
#        frame.GetYaxis().SetLabelSize(0.045)
        frame.GetYaxis().SetTitleSize(0.05)
#        frame.GetYaxis().SetTitleOffset(1.4)
#        frame.GetYaxis().SetTitleFont(42)
#        frame.GetZaxis().SetLabelFont(42)
#        frame.GetZaxis().SetLabelOffset(0.007)
#        frame.GetZaxis().SetLabelSize(0.045)
#        frame.GetZaxis().SetTitleSize(0.05)
#        frame.GetZaxis().SetTitleFont(42)

        if len(units)>0:
            frame.GetXaxis().SetTitle(titlex + " (" +units+")")
            frame.GetYaxis().SetTitle("Events / "+str((hists[0].GetXaxis().GetXmax()-hists[0].GetXaxis().GetXmin())/hists[0].GetNbinsX())+ " "+units)
        else:    
            frame.GetXaxis().SetTitle(titlex)
            frame.GetYaxis().SetTitle("Events")

        frame.Draw()
        stack.Draw("A,HIST,SAME")
        if data !=None:
            dataG.Draw("Psame")              

        ###legend = ROOT.TLegend(0.62,0.6,0.92,0.90,"","brNDC")
        legend = ROOT.TLegend(0.63,0.64,0.93,0.94,"","brNDC") # if using tdrstyle
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)

        legend.SetFillColor(ROOT.kWhite)
        for (histo,label,typeP) in  list(zip(hists,self.labels,self.types))[::-1]:
            if typeP != "data" and typeP !='signal':
                legend.AddEntry(histo,label,"f")
            elif typeP == 'data':
                legend.AddEntry(histo,label,"p")

        for (histo,label,typeP) in  list(zip(hists,self.labels,self.types))[::-1]:
            if typeP == "signal":
                legend.AddEntry(histo,label,"f")

        ###tex_prelim = ROOT.TLatex()
        ###if drawprelim:
        ###    if prelim != "":
        ###        tex_prelim.SetTextSize(0.03)
        ###        tex_prelim.DrawLatexNDC(.11, .91, "#scale[1.5]{CMS}"+" {}".format(prelim))
        ###        tex_prelim.Draw("same")
        
        float_lumi = float(lumi)
        float_lumi = float_lumi/1000.
        ###tex_lumi = ROOT.TLatex()
        ###tex_lumi.SetTextSize(0.035)
        ###tex_lumi.SetTextAlign(31)
        ###tex_lumi.DrawLatexNDC(.93, .91, "13 TeV ({:.1f}".format(float_lumi) + " fb^{-1})")
        ###tex_lumi.Draw("same")
        if outOfFrame:
            CMS_lumi.CMS_lumi(canvas, 4, 0, relPosX=0.077, lumi_13TeV = str(float_lumi), extraText = prelim)
        else:
            CMS_lumi.CMS_lumi(canvas, 4, 10, lumi_13TeV = str(float_lumi), extraText = prelim)

 #       ROOT.SetOwnership(legend,False)

        legend.Draw()
        if self.log:
            canvas.SetLogy()
#       canvas.SetLeftMargin(canvas.GetLeftMargin()*1.15)
        canvas.Update()



        if verbose:
            print("---------------------------")
            print( "Signal = %f" %(signal))
            print( "Bkg    = %f" %(background))
            if data is not None:
                print ("Observed = %f"%(data.Integral()))
                integral = data.IntegralAndError(1,data.GetNbinsX(),error)
                if background>0.0:
                    print ("Data/Bkg= {ratio} +- {err}".format(ratio=integral/background,err=math.sqrt(error.value*error.value/(background*background)+integral*integral*backgroundErr/(background*background*background*background))))

        canvas.RedrawAxis()



        canvas.Update()
        ###plot={'canvas':canvas,'stack':stack,'legend':legend,'data':data,'dataG':dataG,'hists':hists,'prelim':tex_prelim, 'lumi': tex_lumi}
        plot={'canvas':canvas,'stack':stack,'legend':legend,'data':data,'dataG':dataG,'hists':hists}

        return plot


# The nostack option normalizes the background and signal
# contributions separately. Without this all MC contributions
# are normalized together and drawn stacked
    def draw_comp(self,var,cut,model, titlex = "", units = "",expandY=0.0,nostack=True,prelim="Work in progress",SFs = "(1)", outOfFrame = 1): 
        ###canvas = ROOT.TCanvas("canvas","")
        canvas = ROOT.TCanvas("canvas", "", 800, 600)
#        ROOT.gStyle.SetOptStat(0)
#        ROOT.gStyle.SetOptTitle(0)
#        canvas.Range(-68.75,-7.5,856.25,42.5)
#        canvas.SetFillColor(0)
#        canvas.SetBorderMode(0)
#        canvas.SetBorderSize(2)
#        canvas.SetTickx(1)
#        canvas.SetTicky(1)
#        canvas.SetLeftMargin(0.15)
###        canvas.SetRightMargin(0.05)
#        canvas.SetTopMargin(0.05)
#        canvas.SetBottomMargin(0.15)
#        canvas.SetFrameFillStyle(0)
#        canvas.SetFrameBorderMode(0)
#        canvas.SetFrameFillStyle(0)
#        canvas.SetFrameBorderMode(0)

        canvas.SetTopMargin(0.055) # testing tdr style
        canvas.SetBottomMargin(0.12) # testing
        canvas.SetLeftMargin(0.12) # testing


        canvas.cd()
        hists=[]
        labels = {}
        stack = ROOT.THStack("stack","")
        
        signal=0
        background=0
        backgroundErr=0
        
        data=[]

        cutL="("+self.defaultCut+")*("+cut+")"
        scale = 0.0



        for (plotter,typeP,label,name) in zip(self.plotters,self.types,self.labels,self.names):
            hist = plotter.hist1d(var,cutL+"*("+SFs+")","1",model,titlex,units)
            hist.SetFillStyle(0)
            hist.SetName(name+label)
            labels[hist.GetValue()] = label
            if nostack:
                if hist.Integral() == 0:
                    stack.Add(hist.GetValue())
                    hists.append(hist)
                    continue
                hist.Scale(1.0/hist.Integral())
                stack.Add(hist.GetValue())
                hists.append(hist.GetValue())
            else:
                if typeP =="data":
                    if hist.Integral() > 0:
                        hist.Scale(1.0/hist.Integral())
                    data.append(hist.GetValue())
                else:
                    scale += hist.Integral()
                    hists.append(hist.GetValue())
                    
        if nostack:
            stack.Draw("hist,nostack")
        else:
            for h in hists:
                h.Scale(1./scale)
                stack.Add(h)
            stack.Draw("hist")
            for h in data:
                h.Draw("hist,same")

        ###canvas.SetLeftMargin(canvas.GetLeftMargin()*1.15)
        stack.SetMinimum(0)
        if len(units):
            stack.GetXaxis().SetTitle(titlex + " ["+units+"]")
        else:
            stack.GetXaxis().SetTitle(titlex)

        stack.GetYaxis().SetTitle("a.u.")
        stack.GetYaxis().SetTitleOffset(0.9)
        stack.GetYaxis().SetTitleSize(0.05)
        stack.GetXaxis().SetTitleSize(0.05)

        #legend = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
        legend = ROOT.TLegend(0.61,0.64,0.91,0.94,"","brNDC") # if using tdrstyle
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)

        legend.SetFillColor(ROOT.kWhite)
        for histo in labels.keys():
            legend.AddEntry(histo, labels[histo], 'lf')
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        ROOT.SetOwnership(legend, False)
        legend.Draw()
        
        #if drawprelim:
        #    tex_prelim = ROOT.TLatex()
        #    if prelim != "":
        #        tex_prelim.SetTextSize(0.03)
        #        tex_prelim.DrawLatexNDC(.11, .91, "#scale[1.5]{CMS}"+" {}".format(prelim))
        #        tex_prelim.Draw("same")
        
        if outOfFrame:
            CMS_lumi.CMS_lumi(canvas, 0, 0, relPosX=0.077, extraText = prelim)
        else:
            CMS_lumi.CMS_lumi(canvas, 0, 10, extraText = prelim)





        canvas.Update()

        return {'canvas': canvas, 'stack': stack, 'legend': legend, 'data': data, 'hists': hists}

