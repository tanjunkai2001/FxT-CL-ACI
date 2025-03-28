function [x,dx] = func_reference(t)
    global a_ref f_ref
    a = a_ref;
    f = f_ref;
    x = a*[2*sin(f*t); cos(2.3*f*t)]+2*a;
    dx = a*f*[2*cos(f*t); -2.3*sin(2.3*f*t)];