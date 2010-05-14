function val = gauss(x,sigma)
    val =  (1/(sigma^2*sqrt(2*pi)) * exp(-((x^2)/(2*sigma^2))));
