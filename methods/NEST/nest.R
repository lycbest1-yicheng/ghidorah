library(R6)

MACRO_CAPVAR = TRUE

NEST = 
  R6Class("NEST",
    public = list(
        hx = NULL,
        hsig = NULL,
        parK = 5, #number of partition for SURE estimation
        trainingSet = NULL,
        includeVar = MACRO_CAPVAR,
    initialize = function(trainingSet){
        required = c("y_scaled", "scale", "sample_size")
        assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled_A, y_scaled_B, scale, sample_size!")
        df = trainingSet %>% (private$enrichTrainingData) 
        total.var = mean(trainingSet$y_scaled^2)  
        noise.var = var(1/trainingSet$sample_size)
        prior.sd = sqrt(total.var/noise.var)
        parOpt = optim(c(prior.sd, sqrt(noise.var)),private$sureObj, df = df, method="L-BFGS-B", lower=c(0,0), upper = c(5,5), control = list(trace =0))
        self$hx = parOpt$par[1]
        self$hsig = parOpt$par[2]
        self$trainingSet = df
    },
    train = function(){invisible(self)},
    predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        assertthat::assert_that(all(!is.null(self$hx), !is.null(self$hsig)), msg="Parameter nu has not been initiated from training data!")
        df = newdata %>% (private$enrichTestingData)
        pred = mapply(private$predictfn, x=df$y_scaled, sigma=df$sigma, MoreArgs=list(xobs=self$trainingSet$y_scaled, sigmaobs=self$trainingSet$sigma, h=c(self$hx, self$hsig)), SIMPLIFY=TRUE)
        if (self$includeVar) {
          return(list(fit = unlist(pred[1,]), varfit = unlist(pred[2,])))
        }
        return(unlist(pred[1,]))
    }
    ),
    private = list(
      training = NULL,
      enrichTrainingData = function(df) {
        df %>% (private$partitionDf) %>% mutate(sigma = sqrt(1/df$sample_size))
      },
      enrichTestingData = function(df) {
        df %>% mutate(sigma = sqrt(1/df$sample_size))
      },
      partitionDf = function(df) {
        sample = runif(nrow(df))
        groupId = findInterval(sample, seq(0,1,1/self$parK))
        df = df %>% mutate(groupId = groupId)
      },
      sureX = function(x, sigma, groupId, xj, sigmaj, groupIdj, h) { 
        x.diff = (xj - x)/(h[1]*sigmaj)
        x.den = exp(-x.diff^2/2)/sigmaj*(groupIdj!=groupId)
        sigma.den = exp(-(sigma-sigmaj)^2/2/h[2]^2)*(groupIdj!=groupId)
        f.hat = sum(x.den*sigma.den)
        f1d.hat = sum(x.den*sigma.den*x.diff/(h[1]*sigmaj))
        f2d.hat = sum(x.den*sigma.den*(x.diff^2-1)/(h[1]*sigmaj)^2)
        sure.x = sigma^2 + sigma^4*(2*f.hat*f2d.hat-f1d.hat^2)/f.hat^2
        if (is.finite(sure.x)) {
          return(sure.x)
        }else{
          return(1e5) #return some value large
        }
      },
      sureObj = function(par, df) {
        sum(mapply(private$sureX, x=df$y_scaled, sigma=df$sigma, groupId=df$groupId, MoreArgs = list(h = par, xj=df$y_scaled, sigmaj=df$sigma, groupIdj=df$groupId)))
      },
      predictfn = function(x, sigma, xobs, sigmaobs, h) {
        x.diff = (xobs - x)/(h[1]*sigmaobs)
        x.den = exp(-x.diff^2/2)/sigmaobs
        sigma.den = exp(-(sigma-sigmaobs)^2/2/h[2]^2)
        f.hat = sum(x.den*sigma.den)
        f1d.hat = sum(x.den*sigma.den*x.diff/(h[1]*sigmaobs))
        f2d.hat = sum(x.den*sigma.den*(x.diff^2-1)/(h[1]*sigmaobs)^2)
        x.mean = x + sigma^2*f1d.hat/f.hat
        x.var = sigma^2*(1+sigma^2*(f2d.hat*f.hat - f1d.hat^2)/f.hat^2)
        if (x.var <=0) x.var = sigma^2
        return(list(mean = x.mean, var = x.var))
      }
    )

  )