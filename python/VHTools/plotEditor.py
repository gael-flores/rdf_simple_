import ROOT
class plotter:

    def getCanvas(self, title, w = 800, h = 600, grid = False, logy = False):
        c = ROOT.TCanvas(title, '', w, h)
        #c.SetLeftMargin(0.11)
        #c.SetRightMargin(0.05)
        if grid:
            c.SetGrid()
        if logy:
            c.SetLogy()
#        c.SetTickx()
#        c.SetTicky()
        return c

    def getCanvas2D(self, title, w = 800, h = 600):
        c = ROOT.TCanvas(title, '', 800, 600)
        #c.Range(-1.373628,-1.338583,113.6813,111.1024);
        #c.SetFillColor(0);
        #c.SetBorderMode(0);
        #c.SetBorderSize(2);
        c.SetRightMargin(0.16);
        c.SetLeftMargin(0.13)
        #c.SetFrameBorderMode(0);
        #c.SetFrameBorderMode(0);
        return c

    def savePlot(self, c, name):
        c.SaveAs(name+".root")
        c.SaveAs(name+".png")
        c.SaveAs(name+".pdf")

    def drawPrelim(self, c, text = "Preliminary", size = 0.03):
        c.cd()
        prelim = ROOT.TLatex()
        prelim.SetTextSize(size)
        prelim.DrawLatexNDC(.13, .95, "#scale[1.5]{CMS}"+" {}".format(text))
        prelim.Draw("same")

    def drawPrelim_ZH(self, c, text = "Work in progress", size = 0.03):
        c.cd()
        prelim = ROOT.TLatex()
        prelim.SetTextSize(size)
        prelim.DrawLatexNDC(.13, .95, "#scale[1.5]{CMS}"+" {}".format(text))
        prelim.Draw("same")

    def drawCOM(self, c, lumi):
        c.cd()
        com = ROOT.TLatex()
        com.SetTextSize(0.035)
        com.SetTextAlign(31)
        com.SetTextFont(52)
        com.DrawLatexNDC(0.91, 0.95, "13 TeV ({:.1f}".format(lumi/1000.) + " fb^{-1})")
        com.Draw("same")

    def getLegend(self, x1, y1, x2, y2, size = 0.04):
        l = ROOT.TLegend(x1, y1, x2, y2)
        l.SetFillStyle(0)
        l.SetBorderSize(0)
        l.SetTextSize(size)
        return l

    def labelGraph(self, g, xlabel, ylabel, x1, x2, y1, y2, title = ''):
        g.SetTitle(title)
        g.GetXaxis().SetTitle(xlabel)
        g.GetYaxis().SetTitle(ylabel)
        g.GetXaxis().SetTitleOffset(1.2)
        g.GetYaxis().SetTitleOffset(1.4)
        g.GetXaxis().SetRangeUser(x1, x2)
        g.GetYaxis().SetRangeUser(y1, y2)

    def formatGraph(self, g, color, marker = 21, size = 1, width = 2):
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(marker)
        g.SetMarkerSize(size)
        g.SetLineWidth(width)