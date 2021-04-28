%% Function to return musical note and pitch deviation in cents

function [r] = Which_note(f)
    n = 4.75+(log(f)-log(440))./log(2);
    octave = floor(n+1/24);
    m = round((n-octave)*12 + 1);
    nn = ["C","C#","D","Eb","E","F","F#","G","G#","A","Bb","B"];
    r = nn(m);
    for i = 1:length(f)
        r(i) = sprintf("%s%d%+.1f¢",r(i),octave(i),100*((n(i)-octave(i))*12-m(i)+1));
    end
end