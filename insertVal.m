function M = insertVal(matrix,x,y,eq,sym)
    if not(has(eq, sym)) 
        matrix(x,y) = 0; 
    else
        % Trick - równania są liniowe, wiec policzenie pochodnej po zmiennej
        % zwroci wspolczynnik przez który mnożona jest zmienna
        matrix(x,y) = diff(eq,sym);
    end
    M = matrix;
end
