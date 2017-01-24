n = 1000;

%v = [1:n];
%P = perms(v);
%a = P(345600,:)

a = rand(1,n)
b = a;
p = [1:n];

for i= 1:n-1
    for j = i+1:n
        if(a(i) < a(j))
            
            val = a(i);
            a(i) = a(j);
            a(j) = val;
            
            ival = p(i);
            p(i) = p(j);
            p(j) = ival;
        end
    end
end

a
p

for i=1:n
    invp(p(i)) = i;
end
    
invp

i = 1;
val1 = b(p(i));
t = zeros(n,1);
fini = false;

while (not(fini))
  
    while(t(p(i)) == 0)
        t(p(i)) = 1;
        val2 = b(i);
        b(i) = val1;
        val1 = val2;
        i = invp(i);
    end
    
    j = 1;
    while((j <= n) && (t(j) == 1))
        j = j + 1;
    end
    
    if(j > n)
        fini = true;
    else
        i = j;
        val1= b(p(i));
    end
end

b
norm(a-b)
  