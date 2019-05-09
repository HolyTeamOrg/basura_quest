# PROPÓSITO

Se trata de diseñar un mecanismo para coger un cubo de basura y volcarlo dentro del camión.

 Esto se hace con dos barras accionadas mediante un pistón hidráulico que están unidas a los extremos del cubo y al camión. 
 
 Se trata de hallar los puntos del camión donde hay que colocar las barras para que el cubo cumpla unos requisitos durante su desplazamiento:
1. Su inclinación está comprendida entre 60º y 90º hasta que pasa la horizontal del camión.
2. Vuelca dentro del camión con un ángulo de -66º y entre x=133, x=400.
3. No choca con el camión.


## CÓMO USAR BASURA.M

1. Ejecutar.
1. Dar un nº entero de iteraciones >=1.
1. Dar un nº de nodos entero >=1. 

**Cuidado, para números muy grandes el tiempo de análisis se dispara.**


## FUNCIONAMIENTO

El método iterativo usado es similar al método de elementos finitos. El factor de calidad resultante consta de 3 argumentos, que son `(A, B, C)`, donde:
* A: porcentaje del tiempo que cumple el requisito de diseño 1.
* B: lo cerca que está de cumplir el requisito de diseño 2.
* C: Validez total de la mejor solución analizada, entre 0 y 1.

## SOBRE LAS VERSIONES DE BASURA.M

- basura1: Versión inicial, apoyos en el cubo fijos en su lado derecho.
- basura2: Apoyos fijos en el cubo variables en el eje x.
- basura3: Apoyos fijos en el cubo variables en el eje x y su altura en el eje y.
- basura4: Apoyos en el cubo en cualquier punto.