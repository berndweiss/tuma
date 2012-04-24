##' tuma: A graphical user interface for teaching univariate meta-analysis
##'
##' This package provides a graphical user interface (GUI) to visualize
##' univariate fixed effect (FEM) and random effects model (REM) meta-analysis.
##' By interactively manipulating effect sizes and standard errors, it is
##' possible to better understand their relationship with summary (FEM and REM)
##' and heterogeneity related statistics (Q, I^2 and tau^2). Additionally, two
##' plots are created which show a typical forest plot as well as a density plot
##' with N(hat{theta}, tau^2). These plots also change if the individual
##' effects sizes or standard errors are manipulated. This work was supported
##' by a fellowship within the "Postdoc-Programme" of the German Academic
##' Exchange Service (DAAD).
##'
##' The tuma GUI can be started by typing \code{tumaGUI()}.
##'
##' @docType package
##' @name tuma
##' @aliases tuma package-tuma
NULL



.onLoad <- function(...){
    txt <- paste("\nThe GUI can be loaded with the function \"tumaGUI()\".\n",
                 "For more information (e.g. how to use another data set), see ?tumaGUI", sep = "")
    packageStartupMessage(txt)
    initializeData()
    tumaGUI()
}




##' Initialize some sample data sets
##'
##' @keywords internal
initializeData <- function(){

    ## Source: Raudenbush, S. W. (1984). Magnitude of Teacher Expectancy Effects
    ## on Pupil IQ as a Function of the Credibility of Expectancy Induction: A
    ## Synthesis of Findings From 18 Experiments. Journal of Educational
    ## Psychology, 76(1), 85--97.
    raudenbush_tuma <<- data.frame(T = c(0.03, 0.12, -0.14, 1.18, 0.26, -0.06, -0.02,
                                  -0.32, 0.27, 0.80, 0.54, 0.18, -0.02, 0.23, -0.18,
                                  -0.06, 0.30, 0.07, -0.07),
                                  se = c(0.125, 0.147, 0.167, 0.373,  0.369, 0.103,
                                  0.103, 0.220, 0.164, 0.251, 0.302, 0.223, 0.289,
                                  0.290, 0.159, 0.167, 0.139, 0.094,  0.174),
                                  studlab = 1:19)


    ## Source: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
    ## (2009). Introduction to meta-analysis. John Wiley and Sons.
    ## (Table 14.7, p. 98)
    borenstein_tuma <<- data.frame(T = c(0.5493, 0.6931, 0.4236,0.2027,0.8673,0.4847),
                                  se = sqrt(c(0.027,0.0115,0.0455,0.0025,0.0175,0.0213)),
                                  studlab = c("Fonda","Newman","Grant","Granger",
                                  "Milland","Finch"))
}




##' Formats numeric tkentries as x.xxx
##'
##' @param x numeric
##' @return formated string
##' @keywords internal
printStats <- function(x){
    x <- formatC(round(x, digits=3), format="f", digits=3)
    return(x)
}




##' Forces effect size (T) elements in a list of tcl variables to be equal
##'
##' This function is called by pressing the "Equalize" sliders and overwrites
##' individual effect size values in a list of tcl variables.
##' @keywords internal
equalizeT <- function(...){
    k <- length(lstTclVarT)
    ## Check if this is the very first call by tkgrid(...), then do nothing,
    ## otherwise it stops with "Error in `*tmp*`[[i]] : subscript out of bounds"
    if(tclvalue(TclVarEqualT) == "1.000"){
        ## do nothing...
    }else{
        for(i in 1:k){
            tclvalue(lstTclVarT[[i]]) <- as.numeric(tclvalue(TclVarEqualT))
        }
        OnOK()
    }
}



##' Forces standard errors (se) elements in a list of tcl variables to be equal
##'
##' This function is called by pressing the "Equalize" sliders and overwrites
##' individual standard error values in a list of tcl variables.
##' @keywords internal
equalizeSe <- function(...){
    k <- length(lstTclVarT)
    ## Check if this is the very first call by tkgrid(...), then do nothing,
    ## otherwise it stops with "Error in `*tmp*`[[i]] : subscript out of bounds"
    if(tclvalue(TclVarEqualSe) == "0.001"){
        ## do nothing
    }else{
        for(i in 1:k){
            tclvalue(lstTclVarSe[[i]]) <- as.numeric(tclvalue(TclVarEqualSe))
        }
        OnOK()
    }
}




##' Creates a simple forest plot without
##'
##' @param T
##' @param se
##' @param index
##' @param ll
##' @param ul
##' @param mylinewidth
##' @keywords internal
plotForest <- function(T, se, index, ll, ul, mylinewidth){
    w <- 1/(se^2)
    point.size <- sqrt(w/min(w))
    points(T, index+3, pch = 15, cex=point.size)
    segments(ll,index+3,ul,index+3, lwd=mylinewidth)
}




##' Convert list of tclVar variables to vector
##'
##' Takes a list of tclVar objects and returns a numeric vector or a character
##' vector which then can be processed by other functions.
##'
##' @param x a list of tclVar() objects
##' @keywords internal
tclVarList2vector <- function(x){
    xvec <- unlist(lapply(x, function(x)tclvalue(x)))

    ## test if xvec contains any letters from a:z or A:Z
    ## if yes, then xvec contains study labels
    xvec.test <- grep("[a:zA:Z]", xvec)
    if(length(xvec.test) == 0){
        xvec <- as.numeric(xvec)
    }
    return(xvec)
}




##' Open new plot window and creates forest and density plot
##'
##' @param T numeric vector of effect sizes
##' @param se numeric vector of standard errors
##' @param k length of T
##' @param studlab character vector of study labels
##' @param myfontsize font size (cex)
##' @param mylinewidth line width (lwd)
##' @keywords internal
plotMeta <- function(T, se, k, studlab, myfontsize, mylinewidth){

    ## Run univariate meta-analysis
    resMeta <- metagen(TE = T, seTE = se)

    ## Calculate confidence intervals
    ll <- T - 1.96*se
    ul <- T + 1.96*se

    ## Define grafical parameters
    x.min <- round(min(ll)-1)
    x.max <- round(max(ul)+1)
    x.seq <- seq(x.min, x.max, 0.01)

    ## Define layout for plots: 2 rows, 1 column
    layout(matrix(c(1, 2), 2, 1))#, width = c(0.32, 0.68))


    ## Plot 1: Forest plot
    par(mar = c(2, 5, 2, 0.5), cex.axis=myfontsize, cex.lab=myfontsize, yaxt="n")

    plot(NULL, NULL, xlim = c(x.min,x.max), ylim = c(1,k+3),
         main = "Forest plot (FEM)", xlab = "Effect size", ylab = "")

    axis(2, at=1:k+3, labels = FALSE, cex.lab=myfontsize, cex.axis=myfontsize)

    text(y = seq(1, k+3, 1), par("usr")[1],
         labels = c("FEM","REM","",rev(studlab)),
         pos = 2, xpd = TRUE, cex=myfontsize)

    plotForest(rev(T), rev(se), 1:k, rev(ll), rev(ul), mylinewidth)
    abline(v=resMeta$TE.random, col = "red",  lwd = mylinewidth)
    abline(v=resMeta$TE.fixed, col = "black", lty=2,  lwd = mylinewidth)
    abline(v=0, col = "black", lty=1, lwd=mylinewidth)


    ## Plot overall effect size
    plotSummES <- function(ll, x, ul, y, col, lty){
        len <- ul - ll
        hig <- 0.4 * len
        polygon(x=c(ll, ll+0.5*len, ul, ll+0.5*len, ll),
                y=c(y, y+hig, y, y-hig, y), col = col, lty =lty, lwd=mylinewidth)
    }

    plotSummES(resMeta$TE.fixed - 1.96*resMeta$seTE.fixed,
            resMeta$TE.fixed,
            resMeta$TE.fixed + 1.96*resMeta$seTE.fixed,
            1, col = "black", lty=2)

    plotSummES(resMeta$TE.random - 1.96*resMeta$seTE.random,
            resMeta$TE.random,
            resMeta$TE.random + 1.96*resMeta$seTE.random,
            2, col = "red", lty=1)


    ## Plot 2: Densityplot
    par(mar = c(2,5,2,0.5), yaxt="s")

    ## Estimate density and check for invalid (Inf) values
    dens <- dnorm(x.seq, mean = resMeta$TE.random, sd = resMeta$tau)
    if(is.infinite(max(dens))){dens[is.infinite(dens)] <- 0}

    plot(x.seq, dens, ylim=c(0, max(dens)),
         type = "l", xlab = "Effect size", ylab = "Density",
         main = "Distribution of true effect sizes", lwd=mylinewidth)

    points(T,rep(0, k), col="red", cex=myfontsize, pch=16)

    abline(v=resMeta$TE.random, col = "red",  lwd = mylinewidth)
    abline(v=resMeta$TE.fixed, col = "black", lty=2,  lwd = mylinewidth)
    abline(v=0, col = "black", lty=1,  lwd = mylinewidth)

    legend(x=min(x.seq)-min(x.seq)*0.05, bty="n", cex=myfontsize,
           y=max(dens)-(max(dens)*0.05), legend=c("REM","FEM"), lty=c(1,2),
           col=c("red","black"), lwd=mylinewidth)
}




##' Main function to be called by most Tk buttons and sliders
##'
##' Each time a button is pressed, a slider is activated or a new data set is
##' loaded, OnOK() is called and causes the graphics window as well as the summary
##' statistics to be updated.
##'
##' @keywords internal
OnOK <- function(...){

    ## Read fontsize (cex) and linewidt (lwd)
    myfontsize <- as.numeric(tclvalue(fontsize))
    mylinewidth <- as.numeric(tclvalue(linewidth))


    ## Initialize and read T, se and studlab vectors
    k <- length(lstTclVarT)
    T <- tclVarList2vector(lstTclVarT)
    se <- tclVarList2vector(lstTclVarSe)
    studlab <- tclVarList2vector(lstTclVarStudlab)


    ## Creates forest plot and density plot
    plotMeta(T = T, se = se, k = k, studlab = studlab,
             myfontsize = myfontsize, mylinewidth = mylinewidth)


    ## Modify tcl variables for summary output
    resMeta <- metagen(TE = T, seTE = se)
    resMeta <- summary(resMeta)

    tclvalue(tclVarFEM) <- paste("FEM = ", printStats(resMeta$fixed$TE),
                                 "[", printStats(resMeta$fixed$lower), "; ",
                                 printStats(resMeta$fixed$upper), "]")

    tclvalue(tclVarREM) <- paste("REM = ", printStats(resMeta$random$TE),
                                 "[", printStats(resMeta$random$lower), "; ",
                                 printStats(resMeta$random$upper), "]")

    tclvalue(tclVarQ) <- paste("Q = ", printStats(resMeta$Q))

    tclvalue(tclVarQ.df) <- paste("(df =", printStats(resMeta$k-1),
                                  ", p =", printStats(1 - pchisq(resMeta$Q, resMeta$k-1)), ")",
                                  sep = "")

    tclvalue(tclVarI2) <- paste("I2 = ", printStats(resMeta$I2$TE))

    tclvalue(tclVarTau2) <- paste("Tau2 = ", printStats(resMeta$tau^2))
}




##' Print metagen() to standard output
##'
##' @keywords internal
onPrintMAResults <- function(){

    T <- tclVarList2vector(lstTclVarT)
    se <- tclVarList2vector(lstTclVarT)

    resMeta <- metagen(TE = T, seTE = se)
    print(resMeta)
    cat("\n\n")
}




## Initializes and updates the Tcl summary variables (FEM, REM, Q etc. )
##'
##' @keywords internal
updateSummaryStats <- function(){

    T <- tclVarList2vector(lstTclVarT)
    se <- tclVarList2vector(lstTclVarSe)
    k <- length(T)
    resMeta <- summary(metagen(TE = T, seTE = se))

    ##
    tclVarFEM <<- tclVar(paste("FEM = ", printStats(resMeta$fixed$TE),
                               "[", printStats(resMeta$fixed$lower), "; ",
                               printStats(resMeta$fixed$upper), "]"))

    tclVarREM <<- tclVar(paste("REM = ", printStats(resMeta$random$TE),
                               "[", printStats(resMeta$random$lower), "; ",
                               printStats(resMeta$random$upper), "]"))

    tclVarQ <<- tclVar(paste("Q = ", printStats(resMeta$Q)))

    tclVarQ.df <<- tclVar(paste("(df =", printStats(resMeta$k-1),
                               ", p =", printStats(1 - pchisq(resMeta$Q, resMeta$k-1)), ")",
                                sep = ""))

    tclVarI2 <<- tclVar(paste("I^2 = ", printStats(resMeta$I2$TE)))

    tclVarTau2 <<- tclVar(paste("Tau^2 = ", printStats(resMeta$tau^2)))
    }






## -------------------------------------------------------------------------- ##
## +++ GUI related functions                                                  ##
## -------------------------------------------------------------------------- ##


##' Initialize lists of tcl variables containing effect sizes, standard errors
##' and study labels
##'
##' This function initializes several lists of tcl variables. Note that the
##' assignment has to take place in the global environment (i.e. use "<<-").
##' In particular, the following lists of tcl variables are generated:
##' \itemize{
##' \item lstTclVarT
##' \item lstTclVarSe
##' \item lstTclVarStudlab
##' \item lstTclEnT
##' \item lstTclEnSe
##' \item lstTclEnStudlab
##' \item lstTclSlT
##' \item lstTclSlSe
##' }
##'
##' @param T numeric vector of effect sizes
##' @param se numeric vector of standard errors
##' @param studlab character vector of study labels
##' @keywords internal
initTclVarLists <- function(T, se, studlab){
    k <- length(T)

    ## ## define variables of class list
    lstTclVarT <<- vector("list", k)
    lstTclVarSe <<- vector("list", k)
    lstTclVarStudlab <<- vector("list", k)

    lstTclEnT <<- vector("list", k)
    lstTclEnSe <<- vector("list", k)
    lstTclEnStudlab <<- vector("list", k)
    lstTclSlT <<- vector("list", k)
    lstTclSlSe <<- vector("list", k)


    lstTclVarT <<- lapply(T, function(x){tmp <- tclVar(printStats(x))})
    lstTclVarSe <<- lapply(se, function(x){tmp <- tclVar(printStats(x))})

    ## tclVar cannot deal with factors, therefor as.character(studlab)
    lstTclVarStudlab <<- lapply(as.character(studlab),
                                function(x){tmp <- tclVar(x)})
}




##' Creates the Tk graphical user interface
##'
##' drawGUI() creates the graphical user interface. Be aware that all tcl
##' related variables need to be globally assigned via "<<-".
##' @keywords internal
drawGui <- function(){
    tt <<- tktoplevel()

    tktitle(tt) <<- "TUMA: Teaching Univariate Meta-Analysis"

    mytkfontsize <<- 10

    tkfont <<- tkfont.create(size=mytkfontsize)

    fontsize <<- tclVar("1")

    linewidth <<- tclVar("1")

    TclVarEqualT <<- tclVar("1.000")

    TclVarEqualSe <<- tclVar("0.001")


    T <- tclVarList2vector(lstTclVarT)
    se <- tclVarList2vector(lstTclVarSe)
    k <- length(T)
    ## use summary(metagen(...)) to calculate and save univ meta-analysis results
    resMeta <- summary(metagen(TE = T, seTE = se))


    ## Combo box ------------------------------------------------------------ ##
    ## source: http://bioinf.wehi.edu.au/~wettenhall/RTclTkExamples/DropDown.html
    OnOKCB <- function(...){

        ## save value of drop box
        choice <<- meta.datasets[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]

        ## delete all tk objects and plot window
        tkdestroy(tt)
        dev.off()

        ## convert string object to data.frame
        df <- get(choice)
        ## update Tcl variables
        initTclVarLists(T = df$T, se = df$se, studlab = df$studlab)

        ##
        drawGui()
        ##OnOK()
    }

    tmp <- ls(name = globalenv())
    meta.datasets <<- tmp[grep("._tuma", tmp)]
    comboBox <- tkwidget(tt, "ComboBox", editable = FALSE, values = meta.datasets)

    OK.but <- tkbutton(tt, text = " Change data set ", command = OnOKCB)
    tkgrid(comboBox, row = "0", column = "0", columnspan = 4, sticky = "w")
    tkgrid(OK.but, row = "0", column = "3", columnspan = 2, sticky = "we")
    ## Combo box ------------------------------------------------------------ ##


    ## empty row  ----------------------------------------------------------- ##
    tkgrid(tklabel(tt, text = ""))


    ## Print summary stats -------------------------------------------------- ##

    ## Initialize summary stats and heterogeneity tests
    updateSummaryStats()

    ## First row: Overall effect size
    tkgrid(tklabel(tt, text="Overall ES: "), row="2", column="0", columnspan=2,
           sticky = "w")
    ## FEM
    tkLabelFEM <- tklabel(tt, text = tclVarFEM, width = 25)
    tkconfigure(tkLabelFEM, textvariable = tclVarFEM)
    tkgrid(tkLabelFEM, row="2", column="1", columnspan = 2, sticky = "we")
    ## REM
    tkLabelREM <- tklabel(tt, text=tclvalue(tclVarREM), width = 25)
    tkconfigure(tkLabelREM, textvariable = tclVarREM)
    tkgrid(tkLabelREM, row="2", column="3", columnspan = 2, sticky = "w")

    ## Second row: Heterogeneity test
    tkgrid(tklabel(tt, text="Test of hetero."), row="3", column="0", columnspan=3,
           sticky = "w")
    ## Q
    tkLabelQ <- tklabel(tt, text = tclVarQ)
    tkconfigure(tkLabelQ, textvariable = tclVarQ)
    tkgrid(tkLabelQ, row="3", column="2", columnspan=2, sticky = "w")
    ## df and p
    tkLabelQ.df <- tklabel(tt, text=tclvalue(tclVarQ.df))
    tkconfigure(tkLabelQ.df, textvariable = tclVarQ.df)
    tkgrid(tkLabelQ.df, row="3", column="3", columnspan=2, sticky = "w")

    ## Row three: Quantification of heterogeneity
    tkgrid(tklabel(tt, text="Quant of hetero"), row="4", column="0", columnspan=3,
           sticky = "w")
    ## I^2
    tkLabelI2 <- tklabel(tt, text = tclVarI2)
    tkconfigure(tkLabelI2, textvariable = tclVarI2)
    tkgrid(tkLabelI2, row="4", column="2", columnspan=2, sticky = "w")
    ## tau^2
    tkLabelTau2 <- tklabel(tt, text=tclvalue(tclVarTau2))
    tkconfigure(tkLabelTau2, textvariable = tclVarTau2)
    tkgrid(tkLabelTau2, row="4", column="3", columnspan=2, sticky = "w")
    ## Print summary stats -------------------------------------------------- ##


    ## empty row ------------------------------------------------------------ ##
    tkgrid(tklabel(tt, text = ""))


    ## Grid of study labels, effect sizes and standard errors --------------- ##
    colhead1 <- tklabel(tt, text=c("Effect size"), font=tkfont)
    colhead2 <- tklabel(tt, text=c("Standard error"), font=tkfont)
    tkgrid(colhead1, row = "12", columnspan=3, column="1", sticky = "w")
    tkgrid(colhead2, row = "12", columnspan=2, column="3", sticky = "w")

    ## Initialize grid of study labels, effect sizes and standard errors...
    ## (note that the initialization of lstTclEnStudlab etc. is done by
    ## initTclVarLists() )
    for(i in 1:k){
        ## study labels
        lstTclEnStudlab[[i]] <- tkentry(tt,width="15",
                                        textvariable=lstTclVarStudlab[[i]], font=tkfont)

        ## effect sizes
        lstTclEnT[[i]] <- tkentry(tt,width="7",textvariable=lstTclVarT[[i]], font=tkfont)
        lstTclSlT[[i]] <- tkscale(tt,from=-2, to=2,resolution = 0.001,
                                  command = OnOK, orient="horizontal",
                                  variable=lstTclVarT[[i]], showvalue=FALSE)

        ## standard errors
        lstTclEnSe[[i]] <- tkentry(tt,width="7",textvariable=lstTclVarSe[[i]], font=tkfont)
        lstTclSlSe[[i]] <- tkscale(tt, from = 0.001, to = 1, resolution = 0.0001,
                                   command = OnOK, orient="horizontal",
                                   variable=lstTclVarSe[[i]], showvalue=FALSE)
    }
    ## ... and create grid
    for(i in 1:k){
        tkgrid(lstTclEnStudlab[[i]],
               lstTclEnT[[i]],
               lstTclSlT[[i]],
               lstTclEnSe[[i]],
               lstTclSlSe[[i]])
    }
    ## Grid of study labels, effect sizes and standard errors --------------- ##


    ## empty row  ----------------------------------------------------------- ##
    tkgrid(tklabel(tt, text = ""))


    ## Equalize row --------------------------------------------------------- ##
    TclEnEqualT <- tkentry(tt,width="7",textvariable=TclVarEqualT, font=tkfont)
    TclSlEqualT <- tkscale(tt,from=-2, to=2,resolution = 0.001,
                           command=equalizeT, orient="horizontal",
                           variable=TclVarEqualT, showvalue=F)

    TclEnEqualSe <- tkentry(tt,width="7", textvariable = TclVarEqualSe,
                            font = tkfont)
    TclSlEqualSe <- tkscale(tt,from = 0.001, to = 1, resolution = 0.0001,
                            command=equalizeSe, orient="horizontal",
                            variable=TclVarEqualSe, showvalue=F)

    tkgrid(tklabel(tt, text="Equalize:"),TclEnEqualT, TclSlEqualT,
           TclEnEqualSe, TclSlEqualSe)
    ## Equalize row --------------------------------------------------------- ##


    ## empty row ------------------------------------------------------------ ##
    tkgrid(tklabel(tt, text = ""))


    ## Buttons row ---------------------------------------------------------- ##
    entryFontsize <- tkentry(tt, width="3", textvariable=fontsize, font=tkfont)
    entryLinewidth <- tkentry(tt, width="3", textvariable=linewidth, font=tkfont)

    buttonUpdate <- tkbutton(tt, text=" Update ",command=OnOK, font=tkfont)

    buttonReset <- tkbutton(tt, text = "Reset", command = function(){
        initData()
        createTclVars()
        OnOK()
    }, font = tkfont)

    buttonExit <- tkbutton(tt, text=" Exit ",command=function(){
        ## Check if graphic device is open and close it
        if(!is.null(dev.list())){
            dev.off()
        }
        ## Close tk window
        tkdestroy(tt)
    }, font=tkfont)

    buttonPrintMA <- tkbutton(tt, text= " Print " ,command = onPrintMAResults,
                              font = tkfont)

    tkgrid(tklabel(tt, text = ""), tklabel(tt, text = ""),
           tklabel(tt, text = ""),
           tklabel(tt, text = "cex:", font=tkfont),
           tklabel(tt, text = "lwd:", font=tkfont))

    tkgrid(buttonUpdate, buttonExit,buttonPrintMA, entryFontsize, entryLinewidth)
    ## Buttons row ---------------------------------------------------------- ##


    tkfocus(tt)
}


    raudenbush_tuma <<- data.frame(T = c(0.03, 0.12, -0.14, 1.18, 0.26, -0.06, -0.02,
                                  -0.32, 0.27, 0.80, 0.54, 0.18, -0.02, 0.23, -0.18,
                                  -0.06, 0.30, 0.07, -0.07),
                                  se = c(0.125, 0.147, 0.167, 0.373,  0.369, 0.103,
                                  0.103, 0.220, 0.164, 0.251, 0.302, 0.223, 0.289,
                                  0.290, 0.159, 0.167, 0.139, 0.094,  0.174),
                                  studlab = 1:19)





##' Start the Tuma GUI (tuma: A graphical user interface for teaching univariate
##' meta-analysis)
##'
##' The Tuma GUI can be started by typing \code{tumaGUI()}. In its current
##' version, the package comes with two data sets (data.frame objects):
##' \itemize{
##' \item \code{borenstein_tuma} (Source: Borenstein, M., Hedges, L. V., Higgins,
##' J. P. T., & Rothstein, H. R. (2009). Introduction to meta-analysis. John
##' Wiley and Sons. (Table 14.7, p. 98))
##' \item \code{raudenbush_tuma} (Source: Raudenbush, S. W. (1984). Magnitude of
##' Teacher Expectancy Effects on Pupil IQ as a Function of the Credibility of
##' Expectancy Induction: A Synthesis of Findings From 18 Experiments. Journal
##' of Educational Psychology, 76(1), 85--97.)
##' }
##' @param data an optional data frame. If a data frame is provided, it has
##' to have at least the following three elements: T, se, studlab. If the name
##' of a newly created data frame ends with "_tuma" (e.g., "borenstein_tuma",
##' "df_tuma", "mydata_tuma" etc.), it will show up in the drop down menu.
##'
##' @author Bernd Weiss \email{bernd.weiss@@uni-koeln.de}
##' @examples
##' tumaGUI()
##' tumaGUI(raudenbush_tuma)
##' @export
##' @import tcltk meta
tumaGUI <- function(data = borenstein_tuma){
    tclRequire("BWidget")

    initTclVarLists(T = data$T,
                    se = data$se,
                    studlab = as.character(data$studlab))
    drawGui()
}


