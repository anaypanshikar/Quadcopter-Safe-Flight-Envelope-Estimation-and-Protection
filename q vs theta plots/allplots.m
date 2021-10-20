b1 = boundary(b1loss(:,3)*180/pi,b1loss(:,4)*180/pi,0.001);
b2 = boundary(b2loss(:,3)*180/pi,b2loss(:,4)*180/pi,0.001);
f1 = boundary(f1loss(:,3)*180/pi,f1loss(:,4)*180/pi,0.001);
f2 = boundary(f2loss(:,3)*180/pi,f2loss(:,4)*180/pi,0.001);
no = boundary(noloss(:,3)*180/pi,noloss(:,4)*180/pi,0.001);

figure()
p1 = patch(b1loss(b1,3)*180/pi,b1loss(b1,4)*180/pi,"w","EdgeColor","b","LineWidth",0.001);
hold on
p2 = patch(b2loss(b2,3)*180/pi,b2loss(b2,4)*180/pi,"w","EdgeColor","r","LineWidth",0.001);
p3 = patch(f1loss(f1,3)*180/pi,f1loss(f1,4)*180/pi,"w","EdgeColor","g","LineWidth",0.001);
p4 = patch(f2loss(f2,3)*180/pi,f2loss(f2,4)*180/pi,"w","EdgeColor","black","LineWidth",0.001);
p5 = patch(noloss(no,3)*180/pi,noloss(no,4)*180/pi,[0.5843 0.3557 0.9982],"EdgeColor",[0.5843 0.3557 0.9982],"LineWidth",0.001);
scatter(-13.1,0,100,"red","x","LineWidth",2)
hold off
alpha(p1,0.5)
alpha(p2,0.5)
alpha(p3,0.5)
alpha(p4,0.5)
alpha(p5,0.5)
grid on
xlabel("theta (deg)")
ylabel("q (deg/s)")
legend("Backward 1 loss","Backward 2 loss","Forward 1 loss","Forward 2 loss","No loss","Trim","location","best")