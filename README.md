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
   
<p align="center">
  <img src="https://github.com/user-attachments/assets/49cc6dec-dc0c-4894-9ac4-b23ff1aefdb8" alt="Figura 1: Sistema nervioso autónomo (simpático y parasimpático).">
  <br><em>Figura 1: Sistema nervioso autónomo (simpático y parasimpático).</em>
</p>


### Frecuencia cardiaca
Por un lado, en la actividad simpática, se incrementa la FC al estimular los receptores β-adrenérgicos del corazón, acelerando el nodo sinoauricular. Por otro lado, la actividad parasimpática disminuye la frecuencia cardíaca mediante la acción del nervio vago, que libera acetilcolina y ralentiza la conducción eléctrica en el corazón.

El balance entre estas dos actividades regula dinámicamente la frecuencia cardíaca en reposo y durante situaciones de estrés o ejercicio.

### Variabilidad de la frecuencia cardíaca (HRV)
La HRV se refiere a las fluctuaciones en los intervalos R-R, estos corresponden al tiempo entre dos picos R sucesivos del ECG, y reflejan la modulación del sistema nervioso autónomo sobre el corazón. Una HRV alta generalmente indica buena regulación autonómica y estado de salud, mientras que una HRV baja puede estar asociada a estrés, fatiga o disfunción del sistema nervioso autónomo.

Las bandas de frecuencia más utilizadas en el análisis espectral de la HRV organizadas de en orden ascendente de frecuencias son:

- *ULF (Ultra Low Frequency):* < 0.003 Hz (requiere registros largos)

- *VLF (Very Low Frequency):* 0.003–0.04 Hz

- *LF (Low Frequency):* 0.04–0.15 Hz (actividad simpática y parasimpática)

- *HF (High Frequency):* 0.15–0.4 Hz (predominantemente parasimpática, relacionada con la respiración)

### Transformada Wavelet
Es una herramienta matemática utilizada para analizar señales en el dominio del *tiempo* y la *frecuencia* **simultáneamente**. A diferencia de la transformada de Fourier, que representa una señal solo en términos de frecuencias globales, la transformada wavelet puede capturar frecuencias locales en distintos momentos del tiempo, como la HRV, en el dominio tiempo-frecuencia. Descompone una señal en pequeñas ondas llamadas wavelets (ondículas), que están localizadas tanto en el tiempo como en la frecuencia. Estas ondículas pueden estirarse (para analizar componentes de baja frecuencia) o comprimirse (para componentes de alta frecuencia).

#### Usos en señales biológicas:

- Análisis de HRV

- Procesamiento de señales EEG, EMG y ECG

- Detección de eventos transitorios o cambios rápidos en la señal

### Tipos comunes de wavelets:

<p align="center">
  <img src="https://github.com/user-attachments/assets/8433c33e-95ca-4b80-aff1-ba07685b5b12" alt="Figura 2: Tipos comunes de wavelets." width="50%">
  <br><em>Figura 2: Tipos comunes de wavelets.</em>
</p>

- Daubechies (db4, db6, etc.): Usada para señales discretas, generalmente en ECG por su forma similar al complejo QRS.
  
- Morlet: Usada para señales contínuas, útil para análisis de frecuencia continua (CWT), especialmente en HRV.

- Symlets: Una variante más simétrica de Daubechies.

- Coiflets: Buena localización en el tiempo y frecuencia.

### Diagrama de flujo

<p align="center">
  <img src="https://github.com/user-attachments/assets/b0c6341d-0450-4af8-9ac3-3ec89dcc4f8b" alt="Figura 4: Diagrama en blanco para anotaciones o flujo personalizado." width="50%">
  <br><em>Figura 4: Diagrama en blanco para anotaciones o flujo personalizado.</em>
</p>

