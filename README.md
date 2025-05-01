# HRV con transformada wavelet
## Descripción 
En este laboratorio trabajamos en el análisis de la variabilidad de la frecuencia cardíaca (HRV) utilizando la transformada wavelet. El objetivo principal fue ver cómo iba cambiando la frecuencia cardíaca de la persona a lo largo del tiempo, especialmente cuando se inducían cambios que afectan la actividad simpática y parasimpática del sistema nervioso.

Para ello, adquirimos una señal ECG durante 300 segundos usando un DAQ y programamos todo en Python. Primero filtramos la señal con un filtro IIR, luego identificamos los picos R y calculamos los intervalos R-R. A partir de ahí analizamos la HRV tanto en el dominio del tiempo como en el dominio tiempo-frecuencia. En este último usamos la transformada wavelet para ver cómo se comporta la potencia espectral en diferentes bandas de frecuencia a lo largo del tiempo.

## Tener en cuenta:
1. Se usa el módulo AD8232 para la toma del ECG.
2. Frecuencia de muestreo de 2000 Hz para asegurar.
3. La duración de la adquisición es de 5 minutos.
4. Se deben tomar los tiempos para los cuales se induce el cambio en la frecuencia cardiaca.
5. Se debe instalar las librerias:
   + Nidaqmx
   + Csv.
   + Numpy.
   + Pandas.
   + Matplotlib.
   + Butter.
   + Filtfilt.
   + Pywt

6. Los electrodos deben estar bien colocados como lo indica una guía médica.
7. Evitar movimientos bruscos o cambios de posición durante la grabación de la señal.
8. Seleccionar una función wavelet adecuada para señales biológicas.

## Fundamento teórico
### Actividades del sistema nervioso autónomo
El sistema nervioso autónomo (SNA) regula las funciones involuntarias del cuerpo, como la frecuencia cardíaca, la presión arterial, la respiración, etc... que se divide en dos ramas para mantener el equilibrio fisiológico (homeostasis).
-   #### Actividad simpática:  
    Activa la respuesta de la llamada lucha o huída, aumentando la frecuencia cardíaca, la fuerza de contracción del corazón aumenta y dilata las vías respiratorias.
    
-   #### Actividad parasimpática:  
   Promueve el estado de descanso o digestión, disminuyendo la frecuencia cardíaca y facilitando la recuperación del cuerpo.

![image](https://github.com/user-attachments/assets/7c717fa5-0dcc-42df-8b18-f9e91c0c4f1c)

### Frecuencia cardiaca
Por un lado, en la actividad simpática, se incrementa la FC al estimular los receptores β-adrenérgicos del corazón, acelerando el nodo sinoauricular. Por otro lado, la actividad parasimpática disminuye la frecuencia cardíaca mediante la acción del nervio vago, que libera acetilcolina y ralentiza la conducción eléctrica en el corazón.

El balance entre estas dos actividades regula dinámicamente la frecuencia cardíaca en reposo y durante situaciones de estrés o ejercicio.

### Variabilidad de la frecuencia cardíaca (HRV)
La HRV se refiere a las fluctuaciones en los intervalos R-R (tiempo entre dos picos R sucesivos del ECG), y refleja la modulación del SNA sobre el corazón. Una HRV alta generalmente indica buena regulación autonómica y estado de salud, mientras que una HRV baja puede estar asociada a estrés, fatiga o disfunción del SNA.

Las bandas de frecuencia más utilizadas en el análisis espectral de la HRV son:

ULF (Ultra Low Frequency): < 0.003 Hz (requiere registros largos)

VLF (Very Low Frequency): 0.003–0.04 Hz

LF (Low Frequency): 0.04–0.15 Hz (actividad simpática y parasimpática)

HF (High Frequency): 0.15–0.4 Hz (predominantemente parasimpática, relacionada con la respiración)
