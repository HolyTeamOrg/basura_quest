function basura4

% clc
home
clear all
close all


%% Condiciones Iniciales

% CUBO. Dimensiones: ancho=180, alto=190.
silueta_cubo(1,:)=[-140, 40, 40+190i, -140+190i, -140];

% Centros en el cubo.
cubo_c1=95i; % SE ELIGE en coordenadas paramétricas del cubo.
cubo_c2=180+95i;
cubo(1,:)=[cubo_c1-140, cubo_c2-140];

ang_inicial=angle(cubo_c2-cubo_c1);
ang_desvio=pi/2-ang_inicial;
ang_actual=ang_inicial;
ang_final=degtorad(-66); % Ángulo de vuelque

% CAMIÓN
% Geometría del camión
camion=[100, 400, 400+375i, 375i 100];

% Cuadrantes de muestreo iniciales
lado_izq=    0;
lado_der=  200;
lado_inf=  200;
lado_sup=  375;
cuadrante1=[lado_izq, lado_der, lado_inf, lado_sup];
cuadrante2=[lado_izq, lado_der, lado_inf, lado_sup];

% Número de posibilidades analizadas
display('__________________________________________________________________________')
m = input('Introduzca el nº de iteraciones a realizar: ');
n = round(sqrt(input('Aproxime el nº de divisiones por iteración: ')))-1;
if n<=0
    n=1;
end
display(' ');
fprintf('Se usarán %d nodos y se analizarán %d casos por iteración. \nEn total serán consideradas %d posibilidades. \r\n',(n+1)^2,(n+1)^4,(n+1)^4*m);
display('Procesando...');


%% Método de cálculo

tic;

figure;
hold on;
% grid on;
axis equal;
plot([-250 450],[0 0],'k') % Suelo
plot(camion,'b')
plot(silueta_cubo,'r')
plot(cubo,'g')
trapez=1; % Forma trapezoidal del cuadrante

for mallado=1:m
    
    k=1;
    porciento=10;
    if mallado>1
        trapez=0;
    end
    
    % Muestreo puntos de apoyo 1
    for c1y=cuadrante1(3):1/n*(cuadrante1(4)-cuadrante1(3)):cuadrante1(4)
        for c1x=cuadrante1(1)-trapez*100/375*(c1y-375):1/n*(cuadrante1(2)-cuadrante1(1)+trapez*100/375*(c1y-375)):cuadrante1(2)
            
            c1=c1x+i*c1y;
            brazo1=norm(cubo(1,1)-c1);
            ang_brazo1=angle(cubo(1,1)-c1);
            plot(real(c1),imag(c1),'b.')
            
            % Muestreo puntos de apoyo 2
            for c2y=cuadrante2(3):1/n*(cuadrante2(4)-cuadrante2(3)):cuadrante2(4)
                for c2x=cuadrante2(1)-trapez*100/375*(c2y-375):1/n*(cuadrante2(2)-cuadrante2(1)+trapez*100/375*(c2y-375)):cuadrante2(2)
                    
                    inclinacion(1)=1;
                    vuelque=0;
                    
                    c2=c2x+i*c2y;
                    brazo2=norm(cubo(1,2)-c2);
                    plot(real(c2),imag(c2),'b.')
                    
                    iter=1;
                    cubo(iter,2*k-1:2*k)=cubo(1,1:2);
                    theta=0;
                    bloqueo=0;
                    
                    centros(k,:)=[c1 c2];

                    % Condiciones de no invasión del camión, etc                    
                    while ((imag(cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i))<(0-375)/(100-0)*real(cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i))+375 ...
                            || (imag(cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i))>375) ...
                            && imag(cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i))>=0 ...
                            && real(cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i))<400)) ... % Fin restricciones esquina 1.
                            && (imag(cubo(iter,2*k-1)+(180+190i-cubo_c1)*exp((ang_actual-ang_inicial)*i))<(0-375)/(100-0)*real(cubo(iter,2*k-1)+(180+190i-cubo_c1)*exp((ang_actual-ang_inicial)*i))+375 ...
                            || imag(cubo(iter,2*k-1)+(180+190i-cubo_c1)*exp((ang_actual-ang_inicial)*i))>375) ...
                            && real(cubo(iter,2*k-1)+(180+190i-cubo_c1)*exp((ang_actual-ang_inicial)*i))<400 ...
                            && real(cubo(iter,2*k-1)+(190i-cubo_c1)*exp((ang_actual-ang_inicial)*i))<=400 ... 
                            && bloqueo==0 ...
                            && ang_final+ang_desvio>=ang_final
                        
                        % Incremento de ángulo
                        theta=theta-1*(pi/180);
                        iter=iter+1;
                        
                        cubo(iter,2*k-1)=brazo1*cos(ang_brazo1+theta)+real(c1)+i*(brazo1*sin(ang_brazo1+theta)+imag(c1));
                        
                        % Intersecciones del vector cubo con biela 2
                        A = [real(cubo(iter,2*k-1)), imag(cubo(iter,2*k-1))];
                        B = [c2x, c2y];
                        a = brazo2; % Radio del segundo círculo
                        b = norm(cubo_c2-cubo_c1); % Radio del primer círculo
                        c = norm(A-B); % Distancia entre circunferencias
                        cosAlpha = (b^2+c^2-a^2)/(2*b*c);
                        
                        if cosAlpha<=1 && cosAlpha>=-1
                            u_AB = (B - A)/c; % Unit vector from first to second center
                            pu_AB = [u_AB(2), -u_AB(1)]; % Perpendicular vector to unit vector
                            intersect_1 = A + u_AB * (b*cosAlpha) + pu_AB * (b*sqrt(1-cosAlpha^2));
                            intersect_2 = A + u_AB * (b*cosAlpha) - pu_AB * (b*sqrt(1-cosAlpha^2));
                            
                            % Seleccionar la intersección apropiada
                            if norm(intersect_1(1)+intersect_1(2)*i-cubo(iter-1,2*k))<norm(intersect_2(1)+intersect_2(2)*i-cubo(iter-1,2*k))
                                cubo(iter,2*k)=intersect_1(1)+i*intersect_1(2);
                            else
                                cubo(iter,2*k)=intersect_2(1)+i*intersect_2(2);
                            end
                            
                            ang_actual=angle(cubo(iter,2*k)-cubo(iter,2*k-1));
                            
                            silueta_cubo(iter,5*k-4:5*k)=[cubo(iter,2*k-1)-(cubo_c1)*exp((ang_actual-ang_inicial)*i),...
                                cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i),...
                                cubo(iter,2*k-1)+(180+190i-cubo_c1)*exp((ang_actual-ang_inicial)*i),...
                                cubo(iter,2*k-1)+(190i-cubo_c1)*exp((ang_actual-ang_inicial)*i),...
                                cubo(iter,2*k-1)-(cubo_c1)*exp((ang_actual-ang_inicial)*i)];
                            
                            
                            %% Parámetros de calidad
                            % Son dos: hasta subir la altura del camión gira entre 90 y 60 grados;
                            % al final del recorrido vuelca con un ángulo ang_final dado.
                            
                            % Calidad de la inclinación
                            if imag(cubo(iter,2*k-1)+(180-cubo_c1)*exp((ang_actual-ang_inicial)*i))<375
                                if ang_actual+ang_desvio<=pi/2 ...
                                        && ang_actual+ang_desvio>=pi/3
                                    inclinacion(iter)=1;
                                else
                                    inclinacion(iter)=0;
                                end
                            end
                            
                            % Calidad ángulo de vuelque
                            if real(cubo(iter,2*k-1)+(180+190i-cubo_c1)*exp((ang_actual-ang_inicial)*i))>=133 ...
                                    && real(cubo(iter,2*k-1)-(cubo_c1)*exp((ang_actual-ang_inicial)*i))<=400
                                if ang_actual+ang_desvio==ang_final
                                    vuelque=1;
                                else
                                    if 1-abs(ang_actual+ang_desvio-ang_final)/pi>vuelque
                                        vuelque=1-abs(ang_actual+ang_desvio-ang_final)/pi;
                                    end
                                end
                            end
                        else
                            bloqueo=1; % Desechar si se bloquea el mecanismo
                        end
                    end
                    
                    calidad(k,:)=[sum(inclinacion)/length(inclinacion), vuelque, (sum(inclinacion)/length(inclinacion))*vuelque];
                    k=k+1;
                    clear vuelque;
                    clear inclinacion; % Borra para obtener el nuevo tamaño iter en k+1
                    if round(k/((n+1)^4)*100)==porciento
                        fprintf('%3d%% \n',porciento);
                        porciento=porciento+10;
                    end
                    [tam1 tam2]=size(silueta_cubo);
                    if iter<tam1
                        silueta_cubo(iter+1:tam1,5*k-4:5*k)=0; % Borra residuos en silueta_cubo tras varias iteraciones.
                        cubo(iter+1:tam1,2*k-1:2*k)=0;
                    end
                end
            end
        end
    end
    
    % Elijo el mejor caso
    [unused k]=max(calidad(:,3));
    
    % Nuevo mallado con los mejores centros en el medio
    cuadrante1=[real(centros(k,1))-1/n*(cuadrante1(2)-cuadrante1(1)+100/375*(imag(centros(k,1))-375)), ...
        real(centros(k,1))+1/n*(cuadrante1(2)-cuadrante1(1)+100/375*(imag(centros(k,1))-375)), ...
        imag(centros(k,1))-1/n*(cuadrante1(4)-cuadrante1(3)), ...
        imag(centros(k,1))+1/n*(cuadrante1(4)-cuadrante1(3))];
    
    if real(centros(k,1))>=400
        cuadrante1(2)=real(centros(k,1));
    end
    if imag(centros(k,1))>=375
        cuadrante1(4)=imag(centros(k,1));
    end
    
    cuadrante2=[real(centros(k,2))-1/n*(cuadrante2(2)-cuadrante2(1))+100/375*(imag(centros(k,2))-375), ...
        real(centros(k,2))+1/n*(cuadrante2(2)-cuadrante2(1))+100/375*(imag(centros(k,2))-375),...
        imag(centros(k,2))-1/n*(cuadrante2(4)-cuadrante2(3)), ...
        imag(centros(k,2))+1/n*(cuadrante2(4)-cuadrante2(3))];
    
    if real(centros(k,2))>=400
        cuadrante2(2)=real(centros(k,2));
    end
    if imag(centros(k,2))>=375
        cuadrante2(4)=imag(centros(k,2));
    end
    
    fprintf('Iteración %d de %d completada. \r\n',mallado,m);
end


%% Representar el mejor de los casos

[tam1 tam2]=size(silueta_cubo);
plot(real(centros(k,1)),imag(centros(k,1)),'go',real(centros(k,2)),imag(centros(k,2)),'go')
t1=0:0.01:2*pi; % t1 es parámetro ángulo
plot(norm(cubo(1,1)-centros(k,1))*cos(t1)+real(centros(k,1)), norm(cubo(1,1)-centros(k,1))*sin(t1)+imag(centros(k,1)),'k:')
plot(norm(cubo(1,2)-centros(k,2))*cos(t1)+real(centros(k,2)), norm(cubo(1,2)-centros(k,2))*sin(t1)+imag(centros(k,2)),'k:')

title({['Simulación de ',num2str(m),' iteraciones con ',num2str((n+1)^2),' nodos.'],...
    ['Factor de calidad: (',num2str(calidad(k,1),'%.2f'),', ',num2str(calidad(k,2),'%.2f'),', ',num2str(calidad(k,3),'%.2f'),').']})

% plot(real(cubo(iter-1,2*k-1)),imag(cubo(iter-1,2*k-1)),'go')
% plot(cubo(iter-1,2*k),'gd')
for iter=1:tam1
    plot(silueta_cubo(iter,5*k-4:5*k),'r');
    plot(cubo(iter,2*k-1:2*k),'g');
end

fprintf('Los nodos buscados son C1=(%.2f, %.2f) y C2=(%.2f, %.2f). \nEl caso dibujado es %0.0f y su factor de calidad es (%.4f, %.4f, %.4f).',...
    real(centros(k,1)),imag(centros(k,1)),real(centros(k,2)),imag(centros(k,2)),k,calidad(k,1),calidad(k,2),calidad(k,3));


%% Fin

tfinal=toc;
display(' ')
fprintf('Completado en %.0f min %.f s. \n',floor(tfinal/60),rem(tfinal,60));
display('__________________________________________________________________________')
display(' ');

end