function basura

% clc
home
clear all
close all


%% Condiciones Iniciales

% Cubo
% asa=0.72;
asa=1;
cubo(1,:)=[40-180*asa+0i 40-180*asa+190i];
silueta_cubo(1,:)=[cubo(1)-asa*180/190*exp(i*pi/2)*(cubo(2)-cubo(1)),...
    cubo(2)-asa*180/190*exp(i*pi/2)*(cubo(2)-cubo(1)),...
    cubo(2)+(1-asa)*180/190*exp(i*pi/2)*(cubo(2)-cubo(1)),...
    cubo(1)+(1-asa)*180/190*exp(i*pi/2)*(cubo(2)-cubo(1)),...
    cubo(1)-asa*180/190*exp(i*pi/2)*(cubo(2)-cubo(1))];

% Geometría del camión
camion=[100+0i 400+0i 400+375i 0+375i 100+0i];

% Ángulo de vuelque
alfa_f=degtorad(-66);

% Cuadrantes de muestreo iniciales
lado_izquierdo= 0;
lado_derecho=   400;
lado_inferior=  250;
lado_superior=  375;
cuadrante1=[lado_izquierdo, lado_derecho, lado_inferior, lado_superior];
cuadrante2=[lado_izquierdo, lado_derecho, lado_inferior, lado_superior];

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
plot([-250 450],[0 0],'k')
plot(camion,'b')
plot(silueta_cubo,'r')
plot(cubo,'g')
trapez=1;

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
            rho1=norm(cubo(1,1)-c1);
            phase1=angle(cubo(1,1)-c1);
            plot(real(c1),imag(c1),'b.')
            
            % Muestreo puntos de apoyo 2
            for c2y=cuadrante2(3):1/n*(cuadrante2(4)-cuadrante2(3)):cuadrante2(4)
                for c2x=cuadrante2(1)-trapez*100/375*(c2y-375):1/n*(cuadrante2(2)-cuadrante2(1)+trapez*100/375*(c2y-375)):cuadrante2(2)
                    
                    inclinacion(1)=1;
                    vuelque=0;
                    
                    c2=c2x+i*c2y;
                    rho2=norm(cubo(1,2)-c2);
                    plot(real(c2),imag(c2),'b.')
                    
                    iter=1;
                    cubo(iter,2*k-1:2*k)=cubo(1,1:2);
                    theta=0;
                    bloqueo=0;
                    
                    centros(k,:)=[c1 c2];
                    
                    % Condiciones de no invasión del camión, etc
                    while ((imag(cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<(0-375)/(100-0)*real(cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))+375 ...
                            || (imag(cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))>375) ...
                            && imag(cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))>=0 ...
                            && real(cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<400)) ... % Fin restricciones esquina 1.
                            && (imag(cubo(iter,2*k)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<(0-375)/(100-0)*real(cubo(iter,2*k)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))+375 ...
                            || imag(cubo(iter,2*k)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))>375) ...
                            && real(cubo(iter,2*k)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<400 ...
                            && real(cubo(iter,2*k)+(1-asa)*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<=400 ... 
                            && bloqueo==0 ...
                            && angle(cubo(iter,2*k)-cubo(iter,2*k-1))>=alfa_f

                        % Incremento de ángulo
                        theta=theta-1*(pi/180);
                        iter=iter+1;
                        
                        cubo(iter,2*k-1)=rho1*cos(phase1+theta)+real(c1)+i*(rho1*sin(phase1+theta)+imag(c1));
                        
                        % Intersecciones del vector cubo con biela 2
                        A = [real(cubo(iter,2*k-1)) imag(cubo(iter,2*k-1))];
                        B = [c2x c2y];
                        a = rho2; % Radio del segundo círculo
                        b = 190; % Radio del primer círculo
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
                            
                            silueta_cubo(iter,5*k-4:5*k)=[cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)),...
                                cubo(iter,2*k)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)),...
                                cubo(iter,2*k)+(1-asa)*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)),...
                                cubo(iter,2*k-1)+(1-asa)*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)),...
                                cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1))];
                            
                            
                            %% Parámetros de calidad
                            % Son dos: hasta subir la altura del camión gira entre 90 y 60 grados;
                            % al final del recorrido vuelca con un ángulo alfa_f dado.
                            
                            % Calidad de la inclinación
                            if imag(cubo(iter,2*k-1)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<375
                                if angle(cubo(iter,2*k)-cubo(iter,2*k-1))<=pi/2 ...
                                        && angle(cubo(iter,2*k)-cubo(iter,2*k-1))>=pi/3
                                    inclinacion(iter)=1;
                                else
                                    inclinacion(iter)=0;
                                end
                            end
                            
                            % Calidad ángulo de vuelque
                            if real(cubo(iter,2*k)-asa*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))>=133 ...
                                    && real(cubo(iter,2*k-1)+(1-asa)*180/190*exp(i*pi/2)*(cubo(iter,2*k)-cubo(iter,2*k-1)))<=400
                                if angle(cubo(iter,2*k)-cubo(iter,2*k-1))==alfa_f
                                    vuelque=1;
                                else
                                    if 1-abs(angle(cubo(iter,2*k)-cubo(iter,2*k-1))-alfa_f)/pi>vuelque
                                        vuelque=1-abs(angle(cubo(iter,2*k)-cubo(iter,2*k-1))-alfa_f)/pi;
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