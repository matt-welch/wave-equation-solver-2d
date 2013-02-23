int main (int argc, char ** argv[]){

    for (i = 1; i <= n; i++) {
        x = erand48(xi);
        y = erand48(xi);
        if(x*x+y*y <= 1.0) count++;
    }
    pi = 4.0*(double)count/(double)n;
    
    return 0;
}
