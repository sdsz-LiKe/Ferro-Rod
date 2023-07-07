alpha = 0.8;
gamma = 114;

delta = 10^(-4);

syms T;

L_count = 1;

t(1) = 11.01000000000;
w(1) = 10.1;
count = 0;

trig=@(t,t_0,w_0,gamma) (gamma/(2*pi))*(sin(2*pi*t_0)-sin(2*pi*t))+(w_0+gamma*cos(2*pi*t_0))*(t-t_0)-((t-t_0).^2);

quad_left=@(t,t_0,w_0,gamma) (gamma./(2.*pi)).*(sin(2.*pi.*t_0)-1)+(w_0+gamma.*cos(2.*pi.*t_0)).*(t-t_0)-((t-t_0).^2);

quad_right=@(t,t_0,w_0,gamma) (gamma./(2.*pi)).*(sin(2.*pi.*t_0)+1)+(w_0+gamma.*cos(2.*pi.*t_0)).*(t-t_0)-((t-t_0).^2);


for j = 2:1000
    s_1 = vpasolve(quad_left(T,t(j-1),w(j-1),gamma) == 0,T,[-Inf Inf]);
    s_2 = vpasolve(quad_right(T,t(j-1),w(j-1),gamma) == 0,T,[-Inf Inf]);
    s_1 = sort(s_1);
    s_2 = sort(s_2);
    s = [s_1;s_2];
    s = sort(s);
    n = size(s);
    conti = true;
    
    if  t(j-1)<=s(2) && t(j-1)>=s(1)
        for i = double(t(j-1) + 10^(-4): delta :s(2))
            tmin = i;
            tmax = i + delta;
            if (trig(tmin,t(j-1),w(j-1),gamma)<0 && trig(tmax,t(j-1),w(j-1),gamma)>0) || (trig(tmin,t(j-1),w(j-1),gamma)>0 && trig(tmax,t(j-1),w(j-1),gamma)<0) || (trig(tmin,t(j-1),w(j-1),gamma)==0) || (trig(tmax,t(j-1),w(j-1),gamma)==0)
                t(j) = (tmin+tmax)/2;
                w(j) = -alpha*(w(j-1)-2*(t(j)-t(j-1))+gamma*(cos(2*pi*t(j-1))-cos(2*pi*t(j))));
                conti = false;
                
                   if w(j) <= w(j-1)
                       count = count + 1;
                   else 
                       count = 0;
                   end
                   
                   if count >= 51 || w(j) < 0
                        t(j) = fix(t(j-1)) + 1 + (asin(1/(pi*gamma))/(2*pi));
                        w(j) = 0;
                        count = 0;
                        L(L_count) = j;
                        L_count = L_count + 1;
                   end 
                   
                break
            end
        end
    
        if conti
           for i =  double(s(n(1)-1): delta :s(n(1)))
               tmin = i;
               tmax = i + delta;
               if (trig(tmin,t(j-1),w(j-1),gamma)<0 && trig(tmax,t(j-1),w(j-1),gamma)>0) || (trig(tmin,t(j-1),w(j-1),gamma)>0 && trig(tmax,t(j-1),w(j-1),gamma)<0) || (trig(tmin,t(j-1),w(j-1),gamma)== 0) || (trig(tmax,t(j-1),w(j-1),gamma)==0)
                   t(j) = (tmin+tmax)/2;
                   w(j) = -alpha*(w(j-1)-2*(t(j)-t(j-1))+gamma*(cos(2*pi*t(j-1))-cos(2*pi*t(j))));
                                  
                   if w(j) <= w(j-1)
                       count = count + 1;
                   else 
                       count = 0;
                   end
                   
                   if count >= 51 || w(j) < 0
                        t(j) = fix(t(j-1)) + 1 + (asin(1/(pi*gamma))/(2*pi));
                        w(j) = 0;
                        count = 0;
                        L(L_count) = j;
                        L_count = L_count + 1;
                   end
                   
                   conti = false;
                   break
               end
           end
        end
    
    
    else
        for i = double(t(j-1) + 10^(-4): delta :s(n))
            tmin = i;
            tmax = i + delta;
            if (trig(tmin,t(j-1),w(j-1),gamma)<0 && trig(tmax,t(j-1),w(j-1),gamma)>0) || (trig(tmin,t(j-1),w(j-1),gamma)>0 && trig(tmax,t(j-1),w(j-1),gamma)<0) || (trig(tmin,t(j-1),w(j-1),gamma)== 0) || (trig(tmax,t(j-1),w(j-1),gamma)==0)
                t(j) = (tmin+tmax)/2;
                w(j) = -alpha*(w(j-1)-2*(t(j)-t(j-1))+gamma*(cos(2*pi*t(j-1))-cos(2*pi*t(j))));
                conti = false;
                
               if w(j) <= w(j-1)
                   count = count + 1;
               else 
                   count = 0;
               end

               if count >= 51  || w(j) < 0
                    t(j) = fix(t(j-1)) + 1 + (asin(1/(pi*gamma))/(2*pi));
                    w(j) = 0;
                    count = 0;
                    L(L_count) = j;
                    L_count = L_count + 1;
               end
                   
                break
            end
        end 
    end
    
    if conti
        t(j) = fix(t(j-1)) + 1 + (asin(1/(pi*gamma))/(2*pi));
        w(j) = 0;
        L(L_count) = j;
        L_count = L_count + 1;
    end
    disp(j)
    disp(w(j))
end

% x = 1: 0.001:1000;
% y_1 = quad_left(x,t(1),w(1),gamma);
% y_2 = trig(x,t(1),w(1),gamma);
% y_3 = quad_right(x,t(1),w(1),gamma);
% plot(x,y_2)

% function y = trig(t,t_0,w_0,gamma)
%     y = (gamma./(2.*pi)).*(sin(2.*pi.*t_0)-sin(2.*pi.*t))+(w_0+gamma.*cos(2.*pi.*t_0)).*(t-t_0)-((t-t_0).^2);
% end
% 
% function y = quad_left(t,t_0,w_0,gamma)
%     y = (gamma./(2.*pi)).*(sin(2.*pi.*t_0)-1)+(w_0+gamma.*cos(2.*pi.*t_0)).*(t-t_0)-((t-t_0).^2);
% end
% 
% function y = quad_right(t,t_0,w_0,gamma)
%     y = (gamma./(2.*pi)).*(sin(2.*pi.*t_0)+1)+(w_0+gamma.*cos(2.*pi.*t_0)).*(t-t_0)-((t-t_0).^2);
% end
