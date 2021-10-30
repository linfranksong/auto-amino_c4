BEGIN{ isd=0 }
            /DV\/DL, AVERAGES OVER 2000000 STEPS/{ isd=1}
            /NSTEP/{ if(isd == 1 ){nstep = $3; }}
            /EPtot/ {
                if(isd == 1 ) {
                        if( nstep >= nstart ) {eptot = $9;
                               print nstep, eptot;} isd = 0;}}

