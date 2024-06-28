load("trip_dsbn.mat","-mat");
t = array;
disp(t)

load("o_vector.mat","-mat");
o = array;
disp(o)

load("d_vector.mat","-mat");
d = array;
disp(d)

% Change the value of n (no. of bus stops);
n = 231;
b = ones(1,n); 

a = zeros(n,1);
oc = zeros(n,1);
dc = zeros(n,1);
diff = zeros(n,1);
dif = zeros(n,1);
di = zeros(n,1);

again = 1;
     
while (again == 1)
    
    for i = 1:n
        oc(i,1) = 0; 
        for j = 1:n      
            oc(i,1) = oc(i,1) + t(i,j) * b(1,j);
        end 
        oi = o(1,i);
        oci = oc(i,1);
        a(i,1) = double(oi)/double(oci);
        %fprintf('%.2f\t',o(1,i)); 
        %fprintf('%.2f\t',oc(i,1));      
        %fprintf('%.4f\t',a(i,1));     
        %fprintf('\n'); 
    end

    for j = 1:n  
        for i = 1:n       
            t(i,j) = t(i,j)*a(i,1);
        end
    end

    %fprintf('\n');  

    %for i = 1:n 
        %for j = 1:n       
            %fprintf('%.2f\t',t(i,j));         
            %fprintf('\t');
        %end    
        %fprintf('\n');
    %end

    %fprintf('\n');  

    for j = 1:n   
        dc(j,1) = 0;     
        for i = 1:n      
            dc(j,1) = dc(j,1) + t(i,j);
        end  
        b(1,j) = double(d(1,j))/double(dc(j,1));      
        %fprintf('%.2f\t', dc(j));      
        %fprintf('%.4f\t', b(j));     
        %fprintf('\n'); 
    end

    for i=1:n    
        for j=1:n     
            t(i,j) = t(i,j) * b(1,j);          
        end
    end

    %fprintf('\n'); 
    for i = 1:n
        oc(i,1) = 0;     
        for j = 1:n     
            oc(i,1) = oc(i,1) + t(i,j);
        end
        diff(i,1) = o(1,i) - oc(i,1); 
        dif(i,1) = diff(i,1); 
        di(i,1) = abs(dif(i,1));   
        %fprintf('%.2f\t', oc(i));   
        %fprintf('%.4f\t', diff(i));   
        %fprintf('%d', dif(i));   
        %fprintf('%d', di(i));  
        %fprintf('\n');
    end

    if (di(1,1)>1 || di(2,1)>1 || di(3,1)>1 || di(4,1)>1 || di(5,1)>1 || di(6,1)>1|| di(7,1)>1|| di(8,1)>1|| di(9,1)>1|| di(10,1)>1) 
        again = 1;
    else 
        again = 0;       
    end

    fprintf('\n');  

    for i = 1:n    
        for j = 1:n       
            fprintf('%.0f\t',t(i,j));       
            fprintf('\t');
        end     
    fprintf('\n');
    end 

end
    
    
disp(t)

filename = 'final_trip_distribution_mat.csv';
csvwrite(filename, t);
fprintf('Your matrix has been successfully exported as "%s".', filename);




