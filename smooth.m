x = -4:4; y = [0 .15 1.12 2.36 2.36 1.46 .49 .06 0];
cs = spline(x, y);
xx = (-4:0.1:4);
figure
plot(x,y,'o',xx,ppval(cs,xx),'-');