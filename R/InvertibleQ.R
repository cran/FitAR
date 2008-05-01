`InvertibleQ` <-
function(phi){
    all(abs(ARToPacf(phi))<1)
}

