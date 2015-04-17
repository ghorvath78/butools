
function pi = GM1StationaryDistr (B, R, K)

    global BuToolsVerbose;

    pi = GIM1_pi(B,R,'MaxNumComp', K+1, 'Verbose', BuToolsVerbose);
end
