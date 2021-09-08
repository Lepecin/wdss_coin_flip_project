library(shiny)

ui<-fluidPage(
  headerPanel("Coin Flip Randomness Tester"),
  sidebarPanel(
    textInput(inputId = "coinstring","Input Coin Flip Series:"),
    textOutput(outputId = "sampledata"),
    selectInput(inputId = "stringlength","Input Series Length:",c(32,64)),
    sliderInput(inputId = "coinprob","Input Coin Bias:",0,1,0.5,step=NULL,round=FALSE),
    textOutput(outputId = "warning"),
    actionButton(inputId = "comp","Compute p-values")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Main Screen",
        tableOutput(outputId = "datatable"),
        plotOutput(outputId = "KSplot")
        ),
      {tabPanel("User Manual",
        textOutput(outputId = "instruc1"),
        textOutput(outputId = "instruc2"),
        textOutput(outputId = "instruc3"),
        textOutput(outputId = "instruc4"),
        textOutput(outputId = "instruc5"),
        textOutput(outputId = "instruc6"),
        textOutput(outputId = "instruc7"),
        textOutput(outputId = "instruc8")
      )}
    )
  )
)


server<-function(input,output){
  library(pracma)
  library(shiny)
  {
    # Function for converting binary strings into
    # corresponding numeric binary vectors
    bitter<-function(x,s){
      y<-numeric(nchar(x))
      for(i in 1:nchar(x)){
        if(
          (substr(x,i,i)!="1")
          &
          (substr(x,i,i)!="0")
        ){
          y[i]<-0
        }else{
          y[i]<-as.numeric(substr(x,i,i))
        }
      }
      return(y)
    }
    
    
    # Function for finding maximum run length of
    # binary sequence
    argmax<-function(x){
      a<-c(x[1],x)
      b<-c(x,x[length(x)])
      k<-abs(a-b);k[1]<-1;k<-k[c(rep(TRUE,length(k)-1),FALSE)]
      j<-1
      v<-0
      for(i in 1:length(k)){
        if(k[i]==0){
          v<-c(v,0)
        }
        if(k[i]==1){
          while((k[i]-k[i+j]!=0)&(i+j<=length(k))){
            j<-j+1
          }
          v<-c(v,j)
          j<-1
        }
      }
      return(max(v))
    }
    
    # Function for finding number of runs of
    # binary sequence
    argrun<-function(x){
      a<-c(x[1],x)
      b<-c(x,x[length(x)])
      return(sum(abs(a-b))+1)
    }
    
    # Primal probability vector generator
    probvec<-function(n,p){
      a<-numeric(n+1)
      for(i in 1:(n+1)){
        a[i]<-(p^(i-1))*((1-p)^(n+1-i))
      }
      return(a)
    }
    
    # Functions for max run length combinatorics
    # vector generator
    comat<-function(n,h,m,z){
      z[n+1,h+1]<-0
      if((n>=1)&(h>=1)){
        for(i in 0:min(n-1,m)){
          z[n+1,h+1]<-z[n+1,h+1]+z[n-i,h]
        }
      }
      if((n>=m+2)&(h>=m+1)){
        for(i in 1:min(n-1-m,m)){
          z[n+1,h+1]<-z[n+1,h+1]-z[n-m-i,h-m]
        }
      }
      if((h==0)&(n<=m)){
        z[n+1,h+1]<-z[n+1,h+1]+1
      }
      if((h==m+1)&(m+1<=n)&(n<=2*m+1)){
        z[n+1,h+1]<-z[n+1,h+1]-1
      }
      return(z)
    }
    filler<-function(n,m){
      z<-matrix(numeric((n+1)^2),n+1,n+1)
      for(i in 0:n){
        for(j in 0:i){
          z<-comat(i,j,m,z)
        }
      }
      return(z)
    }
    
    # CDF and p-value for maximum run length
    maxruncdf<-function(n,p){
      a<-numeric(n+1)
      for(i in 1:n){
        a[1+i]<-sum(filler(n,i)[n+1,]*probvec(n,p))
      }
      return(a)
    }
    maxrunpval<-function(n,p,k){
      a<-maxruncdf(n,p)
      B<-0
      i<-1
      if(a[k+1]<=0.5){
        A<-a[k+1]
        while(B<=A){
          B<-1-a[n+1-i]
          i<-i+1
        }
        B<-1-a[n+3-i]
      }else{
        if(1-a[k]<=0.5){
          A<-1-a[k]
          while(B<=A){
            B<-a[i+1]
            i<-i+1
          }
          B<-a[i-1]
        }else{
          A<-1
        }
      }
      return(A+B)
    }
    
    # Functions for run number combinatorics
    # vector generator
    comat2<-function(r,h,a,b){
      t<-nrow(a)
      if((r%%2==1)&(h>=1)&(r<=t)){
        b[r,h+1]<-b[r,h+1]+a[r,h]
      }
      if((r%%2==1)&(h>=1)&(r>=2)){
        b[r,h+1]<-b[r,h+1]+a[r-1,h]
      }
      if((r%%2==0)&(h<=t)&(r<=t)){
        b[r,h+1]<-b[r,h+1]+a[r,h+1]
      }
      if((r%%2==0)&(h<=t)&(r>=2)){
        b[r,h+1]<-b[r,h+1]+a[r-1,h+1]
      }
      return(b)
    }
    filler2<-function(a){
      t<-nrow(a)
      b<-matrix(numeric((t+1)*(t+2)),t+1,t+2)
      for(i in 1:(t+1)){
        for(j in 0:(t+1)){
          b<-comat2(i,j,a,b)
        }
      }
      return(b)
    }
    distro<-function(n){
      a<-matrix(c(0,1),1,2)
      if(n>1){
        for(i in 1:(n-1)){
          a<-filler2(a)
        }
      }
      b<-a[, rev(seq_len(ncol(a)))]+a
      return(b)
    }
    
    # CDF and p-value for number of runs
    numruncdf<-function(n,p){
      a<-numeric(n+1)
      b<-distro(n)
      for(i in 1:n){
        a[i+1]<-sum(b[i,]*probvec(n,p))
      }
      return(cumsum(a))
    }
    numrunpval<-function(n,p,k){
      a<-numruncdf(n,p)
      B<-0
      i<-1
      if(a[k+1]<=0.5){
        A<-a[k+1]
        while(B<=A){
          B<-1-a[n+1-i]
          i<-i+1
        }
        B<-1-a[n+3-i]
      }else{
        if(1-a[k]<=0.5){
          A<-1-a[k]
          while(B<=A){
            B<-a[i+1]
            i<-i+1
          }
          B<-a[i-1]
        }else{
          A<-1
        }
      }
      return(A+B)
    }
    
    # P-vector generator
    pvecgen<-function(x,p){
      n<-length(x)
      o<-hadamard(n)%*%x
      z<-abs(
        (o-c(n*p,rep(0,n-1)))
        /
          (sqrt(n*p*(1-p)))
      )
      k<-numeric(n)
      for(i in 1:length(z)){
        k[i]<-2*(1-pnorm(z[i]))
      }
      return(k)
    }
    
    # Measure of uniformity
    uniftest<-function(k){
      v<-seq(1/(length(k)),1,1/(length(k)))
      s<-sort(k)
      return(2*abs(1-(s%*%v)/(s%*%s)))
    }
    
    unifplot<-function(k){
      v<-seq(1/(length(k)),1,1/(length(k)))
      plot(sort(k),v,main="Kolmogorov-Smirnov plot for WH p-vector",xlab="Ordered p-vector"
           ,ylab="Quantiles",asp=1)
      abline(a=0,b=1)
      abline(a=0,b=(sort(k)%*%v)/(sort(k)%*%sort(k)),col="red")
    }
    
    pvaldata<-function(y,p){
      d<-data.frame(
        "Tests" = c(
          "Bernoulli dist. test",
          "Longest run test",
          "General runs test",
          "Walsh-Hadamard test"
        ),
        "Pvalues" = c(
          berpval(y,p),
          maxrunpval(length(y),p,argmax(y)),
          numrunpval(length(y),p,argrun(y)),
          NA
        ),
        "Metrics" = c(
          NA,
          NA,
          NA,
          uniftest(pvecgen(y,p))
        )
      )
      return(d)
    }
    
    # P-value generator for
    # Bernoulli test
    berpval<-function(x,p){
      m<-sum(x)/length(x)
      return(2*(1-pnorm(abs(m-p)/sqrt(p*(1-p)/length(x)))))
    }
  }
  
  output$sampledata<-renderText(
  if(nchar(input$coinstring)!=as.numeric(input$stringlength)){
      paste(
        "Input string length is",
        as.character(nchar(input$coinstring)),
        "instead of",
        as.character(input$stringlength)
      )
    }else{
    y<-bitter(input$coinstring,as.numeric(input$stringlength))
      paste(
        "Estimated probability of head:",
        as.character(sum(y)/length(y)),
        "// Max run length:",
        as.character(argmax(y)),
        "// Total number of runs:",
        as.character(argrun(y))
      )
    })
  output$warning<-renderText(
  if((as.numeric(input$coinprob)==1)|(as.numeric(input$coinprob)==0)){
    paste("Hypothesised coin bias can only be non integer")
  })
  pcalc<-eventReactive(input$comp,{
    y<-bitter(input$coinstring,as.numeric(input$stringlength))
    pvaldata(y,as.numeric(input$coinprob))})
  pplot<-eventReactive(input$comp,{
    y<-bitter(input$coinstring,as.numeric(input$stringlength))
    unifplot(pvecgen(y,as.numeric(input$coinprob)))
   })
  output$datatable<-renderTable({
    if(
      (nchar(input$coinstring)==as.numeric(input$stringlength))&
      (as.numeric(input$coinprob)!=0)&
      (as.numeric(input$coinprob)!=1)
      ){
    pcalc()
    }
  })
  output$KSplot<-renderPlot({
    if(
      (nchar(input$coinstring)==as.numeric(input$stringlength))&
      (as.numeric(input$coinprob)!=0)&
      (as.numeric(input$coinprob)!=1)
    ){
      pplot()
    }
  })
  {
  output$instruc1<-renderText({
    "This app measures the user’s ability to simulate a believable series of coin flips."
  })
  output$instruc2<-renderText({
    "Type a series of ones and zeros into the box saying 'Input Coin Flip Series'."
  })
  output$instruc3<-renderText({
      "The ones and zeros can be interpreted as heads and tails of a hypothetical coin respectively."
  })
  output$instruc4<-renderText({
    "The input series is accepted if its length coincides with the chosen sample size, which can be changed in the box saying 'Input Series Length'."
  })
  output$instruc5<-renderText({
    "Then, choose the hypothesised bias of the coin you simulated such that it lies between zero and one using the 'Input Coin Bias' slider."
  })
  output$instruc6<-renderText({
    "Upon pressing the 'Compute p-values' button, the 'Main Screen' tab will show a table of values and a plot."
  })
  output$instruc7<-renderText({
    "The 'Pvalues' column shows the p-values from each hypothesis test applied to the input series, while the 'Metrics' column shows test values that are not p-values."
  })
  output$instruc8<-renderText({
    "Overall, the closer each p-value is to one an the closer the metric is to zero, the more your input looks like a series of coin flips."
  })
}
}
shinyApp(ui=ui,server=server)
