# HRV con transformada Wavelet
## Descripción 
En este laboratorio trabajamos en el análisis de la variabilidad de la frecuencia cardíaca (HRV) utilizando la transformada Wavelet. El objetivo principal fue ver cómo iba cambiando la frecuencia cardíaca de la persona a lo largo del tiempo, especialmente cuando se inducían cambios que afectan la actividad simpática y parasimpática del sistema nervioso.

Para ello, adquirimos una señal ECG durante 300 segundos usando un DAQ y programamos todo en Python. Primero filtramos la señal con un filtro IIR, luego identificamos los picos R y calculamos los intervalos R-R. A partir de ahí analizamos la HRV tanto en el dominio del tiempo como en el dominio tiempo-frecuencia. En este último usamos la transformada Wavelet para ver cómo se comporta la potencia espectral en diferentes bandas de frecuencia a lo largo del tiempo.

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
8. Seleccionar una función Wavelet adecuada para señales biológicas.

## Fundamento teórico
### Actividades del sistema nervioso autónomo
El sistema nervioso autónomo (SNA) regula las funciones involuntarias del cuerpo, como la frecuencia cardíaca, la presión arterial, la respiración, etc... que se divide en dos ramas para mantener el equilibrio fisiológico (homeostasis) [1].

-   #### Actividad simpática:  
    Activa la respuesta de la llamada lucha o huída, aumentando la frecuencia cardíaca, la fuerza de contracción del corazón aumenta y dilata las vías respiratorias [2]
    
-   #### Actividad parasimpática:  
   Promueve el estado de descanso o digestión, disminuyendo la frecuencia cardíaca y facilitando la recuperación del cuerpo [2]
   
<p align="center">
  <img src="https://github.com/user-attachments/assets/49cc6dec-dc0c-4894-9ac4-b23ff1aefdb8" alt="Figura 1: Sistema nervioso autónomo (simpático y parasimpático).">
  <br><em>Figura 1: Sistema nervioso autónomo (simpático y parasimpático).</em>
</p>


### Frecuencia cardiaca
Por un lado, en la actividad simpática, se incrementa la FC al estimular los receptores β-adrenérgicos del corazón, acelerando el nodo sinoauricular. Por otro lado, la actividad parasimpática disminuye la frecuencia cardíaca mediante la acción del nervio vago, que libera acetilcolina y ralentiza la conducción eléctrica en el corazón.

El balance entre estas dos actividades regula dinámicamente la frecuencia cardíaca en reposo y durante situaciones de estrés o ejercicio.

### Variabilidad de la frecuencia cardíaca (HRV)
La HRV se refiere a las fluctuaciones en los intervalos R-R, estos corresponden al tiempo entre dos picos R sucesivos del ECG, y reflejan la modulación del sistema nervioso autónomo sobre el corazón. Una HRV alta generalmente indica buena regulación autonómica y estado de salud, mientras que una HRV baja puede estar asociada a estrés, fatiga o disfunción del sistema nervioso autónomo [3].

Las bandas de frecuencia más utilizadas en el análisis espectral de la HRV organizadas de en orden ascendente de frecuencias son:

- *ULF (Ultra Low Frequency):* < 0.003 Hz (requiere registros largos)

- *VLF (Very Low Frequency):* 0.003–0.04 Hz

- *LF (Low Frequency):* 0.04–0.15 Hz (actividad simpática y parasimpática)

- *HF (High Frequency):* 0.15–0.4 Hz (predominantemente parasimpática, relacionada con la respiración)

### Transformada Wavelet
Es una herramienta matemática utilizada para analizar señales en el dominio del *tiempo* y la *frecuencia* **simultáneamente**. A diferencia de la transformada de Fourier, que representa una señal solo en términos de frecuencias globales, la transformada Wavelet puede capturar frecuencias locales en distintos momentos del tiempo, como la HRV, en el dominio tiempo-frecuencia. Descompone una señal en pequeñas ondas llamadas Wavelets (ondículas), que están localizadas tanto en el tiempo como en la frecuencia. Estas ondículas pueden estirarse (para analizar componentes de baja frecuencia) o comprimirse (para componentes de alta frecuencia) [4]

#### Usos en señales biológicas:

- Análisis de HRV

- Procesamiento de señales EEG, EMG y ECG

- Detección de eventos transitorios o cambios rápidos en la señal

### Tipos comunes de Wavelets:

<p align="center">
  <img src="https://github.com/user-attachments/assets/8433c33e-95ca-4b80-aff1-ba07685b5b12" alt="Figura 2: Tipos comunes de wavelets." width="50%">
  <br><em>Figura 2: Tipos comunes de Wavelets.</em>
</p>

- *Daubechies:* Usada para señales discretas, generalmente en ECG por su forma similar al complejo QRS.
  
- *Morlet:* Usada para señales contínuas, útil para análisis de frecuencia continua (CWT), especialmente en HRV.

- *Symlets:* Una variante más simétrica de Daubechies.

- *Coiflets:* Buena localización en el tiempo y frecuencia.

### Diagrama de flujo
Planteamos una serie de pasos para llevar a cabo el desarrollo del objetivo, los cuales representamos en el siguiente diagrama de flujo:

<p align="center">
  <img src="https://github.com/user-attachments/assets/b0c6341d-0450-4af8-9ac3-3ec89dcc4f8b" alt="Figura 4: Diagrama en blanco para anotaciones o flujo personalizado." width="50%">
  <br><em>Figura 4: Diagrama correspondiente a los pasos a seguir para desarrollar el laboratorio.</em>
</p>

### Adquisición de la señal ECG
Se configura una tarea para el modulo DAQ para que reciba los datos en una frecuencia de muestreo de 2000 Hz, duración de 300 segundos, correspondientes a los 5 minutos.

```python
import nidaqmx
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
from IPython.display import display, clear_output
import csv

# Parámetros de adquisición
puerto= "Dev2/ai0"  # Cambia esto según el canal de tu DAQ
fmuestreo = 2000  # Frecuencia de muestreo en Hz
t_ventana = 2000  # Cantidad de muestras a mostrar en la gráfica
duracion =  300 # Duración de la adquisición en segundos

# Crear estructura de datos para almacenar la señal en tiempo real
esg_buffer = deque([0] * t_ventana, maxlen=t_ventana)  # Guarda los últimos datos capturados
tiempo_buffer = deque(np.linspace(0, (t_ventana - 1) / fmuestreo,t_ventana), maxlen=t_ventana)  # Eje X desde 0

# Archivo de salida para guardar los datos
filename = "esg_data.csv"
with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Tiempo (s)", "Voltaje (V)"])  # Escribir encabezado

    # Configurar la tarea de adquisición
    with nidaqmx.Task() as task:
        task.ai_channels.add_ai_voltage_chan(puerto)
        task.timing.cfg_samp_clk_timing(fmuestreo)

        print("Adquiriendo datos en tiempo real... Presiona 'Ctrl + C' para detener.")

        # Inicializar la figura
        plt.ion()  # Modo interactivo
        fig, ax = plt.subplots(figsize=(10, 4))
        line, = ax.plot(tiempo_buffer, esg_buffer, label="Señal ECG")
        ax.set_xlabel("Tiempo (s)")
        ax.set_ylabel("Voltaje (V)")
        ax.set_title("Señal ECG en Tiempo Real")
        ax.legend()
        ax.grid(True)
        ax.set_ylim([-5, 5])  # Ajuste de voltaje

        # Bucle de adquisición en tiempo real
        try:
            total_samples = int(duracion * fmuestreo)
            for i in range(0, total_samples, t_ventana):
                # Leer múltiples muestras a la vez
                ecg_data = task.read(number_of_samples_per_channel=t_ventana)

                if isinstance(ecg_data, list):  # Confirmar que es una lista
                    ecg_buffer.extend(ecg_data)

                    # Guardar los datos en el archivo
                    for j in range(len(ecg_data)):
                        tiempo = (i + j) / fmuestreo
                        writer.writerow([tiempo, ecg_data[j]])

                    # Actualizar gráfico en tiempo real
                    line.set_ydata(ecg_buffer)
                    clear_output(wait=True)
                    display(fig)
                    plt.pause(0.01)

        except KeyboardInterrupt:
            print("Adquisición detenida por el usuario.")

        plt.ioff()
        plt.show()

print(f"Los datos se han guardado exitosamente en '{filename}'.")
```
Esta parte del código permite almacenar los datos de tiempo y voltaje correspondientemente por 300 segundos para posteriormente guardarlos en un archivo **.csv** y poderlos graficar.

### Filtrado IIR de la señal
Una vez guardada la señal se procede a aplicar un filtro IIR tipo Butterworth pasa banda de 4° orden entre 0.5 y 40 Hz para la señal ECG.

``` python

# --- Cargar señal ---
ruta_csv = 'ecg_data.csv'
df = pd.read_csv(ruta_csv)

# Seleccionar segunda columna correspondiente al voltaje
ecg_signal = df.iloc[:, 1].values

# --- Parámetros ---
fs = 1000  # Hz
lowcut = 0.5
highcut = 40
order = 4  # Orden del filtro

# --- Filtro digital Butterworth ---
low = lowcut / (fs/2)
high = highcut / (fs/2)

b_d, a_d = butter(order, [low, high], btype='band')

# --- Filtrar la señal ---
ecg_filtered = lfilter(b_d, a_d, ecg_signal) #Filtro IIR

# --- Graficar señal original y señal filtrada---
t = np.arange(len(ecg_signal)) / fs

plt.figure(figsize=(20,6))
plt.subplot(2,1,1)
plt.plot(t, ecg_signal)
plt.title('Señal ECG Original')
plt.grid()
plt.ylabel('Amplitud')

plt.subplot(2,1,2)
plt.plot(t, ecg_filtered, color='red')
plt.title('Señal ECG Filtrada (Pasa Banda 0.5-40 Hz)')
plt.grid()
plt.xlabel('Tiempo [s]')
plt.ylabel('Amplitud')

plt.tight_layout()
plt.show()
```
La función lfilter() aplica ese filtro IIR a la señal ECG. A continuación se muestran las gráficas resultantes para la señal sin filtrar y filtrada, correspondientemente.

<p align="center">
  <img src="https://github.com/user-attachments/assets/d859fd90-0508-4b7c-9e3f-81f9970d6599" alt="Figura 5" width="800">
  <br>
  <em>Figura 5: Señal ECG original y señal filtrada con IIR.</em>
</p>


### Detección de picos R
Ahora con la función *find_peaks* detectamos los picos R con un umbral de prominencia y una distancia mínima entre latidos (250 ms). A partir de esos picos, se calculan los intervalos RR, es decir, el tiempo entre latidos consecutivos.

```python
# --- Detección de picos (ajustada) ---
# Reducimos distancia mínima entre picos para mayor sensibilidad
min_distance = int(0.25 * fs)  # 0.25 segundos (250 ms)
picos, _ = find_peaks(ecg_filtered, distance=min_distance, prominence=0.5)

# --- Calcular intervalos entre picos ---
intervalos = np.diff(picos) / fs  # Diferencias de tiempo en segundos
#print ( intervalos)
# --- Graficar ---
t = np.arange(len(ecg_signal)) / fs
```
A partir de estos cálculos graficamos tanto la señal filtrada en donde se muestran los picos, como los intervalos entre picos:

<p align="center">
  <img src="https://github.com/user-attachments/assets/fcf863c7-590a-4938-b5cf-8e6a5153c358" alt="Figura 6" width="700">
  <br>
  <em>Figura 6: Señal ECG filtrada mostrando los picos R.</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/1efcb229-df8f-4511-aa81-ef1a82b07fcf" alt="Figura 5" width="700">
  <br>
  <em>Figura 5: Señal con intervalos RR.</em>
</p>


Después, analizaremos la HRV, la variabilidad de la frecuencia cardíaca en el dominio del tiempo, calculando métricas como la media RR, SDNN (desviación estándar de los intervalos RR), RMSSD (variabilidad a corto plazo), y pNN50 (porcentaje de intervalos consecutivos). Los resultados se muestran en la siguiente tabla:
| Métrica                        | Valor       | Unidad   |
|-------------------------------|-------------|----------|
| Media de RR (Mean HRV)        | 0.6861      | s        |
| Desviación estándar de RR     | 0.2369      | s        |
| SDNN                          | 0.2369      | s        |
| RMSSD                         | 0.2035      | s        |
| pNN50                         | 72.54       | %        |
| Media de los intervalos RR    | 0.6861      | s        |

### Aplicación Transformada Wavelet
Para realizar la transformada Wavelet continua se debe utilizar la libreria de pywt en python para facilitarnos el analisis de HRV por medio de espectograma el cual nos permite una caracterización especifica y detallada de la actividad cardíaca, tambien utilizamos Wavelet de Morlet porque nos permite ver cómo varían las bandas baja y alta frecuencia a lo largo del tiempo y comparte similitudes en las oscilaciones suaves de la señal biologica, se utilizan 128 escalas por que nos proporciona una resolución adecuada para observar la dinámica temporal de las bandas LF y HF en el análisis de HRV, respetando las limitaciones computacionales y manteniendo la interpretación visual clara y precisa.
```python
# RR como señal base para análisis 
rr_times = np.cumsum(intervalos)  # tiempo relativo de cada latido
rr_signal = intervalos - np.mean(intervalos)  # centrar la señal para análisis

# Wavelet CWT 
# Elección de wavelet biológica y escalas
wavelet = 'cmor1.5-1.0'  # Morlet complejo, buena resolución en tiempo-frecuencia
scales = np.arange(1, 128)  # Escalas 

# Aplicar CWT
coef, freqs = pywt.cwt(rr_signal, scales, wavelet, sampling_period=np.mean(intervalos))
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/77172230-4e8f-42e7-8a47-c112a065d757" alt="Figura 6. Espectrograma transformada Wavelet tipo Morlet" width="600">
</p>

<p align="center"><b>Figura 6.</b> Espectrograma transformada Wavelet tipo Morlet</p>


| Analisis temporal y espectral  | 0-100s     | 100-200s   |200-300s|
|--------------------------------|-------------|----------|----------|
|  Banda de alta frecuencia (0.15–0.4 Hz)  |  En esta zona se puede visualizar una alta *actividad parasimpatica* indicando un estado de relajación especialmente alrededor de los 0.35–0.4 Hz, donde se manifiesta con tonalidades amarillas y verdosas  | Se puede notar que se registra una caída progresiva de la potencia espectral en esta banda. La tonalidad cambia hacia el azul y se puede interpretar como una *disminución de la actividad parasimpática*.| La actividad vagal se ve suprimida y se torna a un estado de estrés o *activación simpática*.|
| Banda de baja frecuencia (0.04–0.15 Hz)  | Se presenta una potencia moderada en esta banda, reflejada por algunos patrones verdosos, aunque no tan intensos como en la banda alta frecuencia.| A medida que avanza el tiempo, se observa un leve incremento en la potencia relativa en la banda de baja frecuencia, especialmente alrededor de los 200 s, lo cual puede ser interpretado como una activación simpática relativa, dado que el componente baja frecuencia puede reflejar tanto simpático como parasimpático.|


#### Diferencias entre el análisis en el dominio del tiempo y el dominio tiempo-frecuencia

El análisis en el dominio del tiempo permite observar el cambio  de los intervalos RR a lo largo del tiempo, pero no da información detallada sobre los componentes de frecuencia que conforman la señal.Aunque se puede ver que la señal cambia con el tiempo, no se puede saber con claridad qué frecuencias están presentes ni en qué momento aparecen o desaparecen las frecuencias de baja y alta frecuencia (LF y HF).
En cambio, con el análisis tiempo-frecuencia ,como el que se hace con la Transformada Wavelet Continua (CWT), se puede observar qué frecuencias hay en cada momento. Esto es muy útil en señales como la HRV, que cambian con el estado del cuerpo. 

#### Efecto del uso de diferentes wavelets
La función wavelet que se utilizó afecta mucho los resultados del análisis. En este caso, se usó la wavelet de Morlet, que es efectiva para el analisis las oscilacione suaves y continuos en frecuencias, justo como los que aparecen en la HRV.

Si se utilizaran otras wavelets, como Daubechies, podríamos detectar cambios más bruscos en la señal, pero se vería con menos claridad las frecuencias, lo cual dificultaria identificar con precisión las bandas espectrales. Otra opción sería Mexican Hat, que resalta los picos o transiciones rápidas, pero no muestra bien las oscilaciones suaves.

Por tanto, según el tipo de señal y lo que se quiere analizar, se debe elegir la wavelet.Para el caso de HRV, donde se quiere observar cómo varían las bandas LF y HF a lo largo del tiempo, la wavelet de Morlet es una de las mejores opciones.

#### Aplicaciones reales de esta práctica

Este tipo de análisis tiene muchas aplicaciones útiles en la vida real. Por ejemplo, se puede usar para estudiar cómo está funcionando el sistema nervioso autónomo, identificando si la persona está en un estado más activo, como cuando está estresada, o más relajado. También es muy útil para monitorear el estrés, la fatiga o la recuperación, algo que se aplica mucho en deportistas o personas que hacen trabajos exigentes. Además, puede ayudar a detectar problemas cardíacos o neurológicos de forma temprana. Hoy en día, esta técnica también se usa en relojes inteligentes y otros dispositivos que miden el nivel de estrés, el sueño o las emociones. Incluso en investigaciones sobre el cerebro y la salud mental, el cambio de la frecuencia cardíaca es una señal importante que puede mostrar cómo responde una persona ante distintas situaciones.

## Bliografia 
[1] Coon, E. (2023, julio 3). Introducción al sistema nervioso autónomo. Manual MSD versión para público general; Manuales MSD. https://www.msdmanuals.com/es/hogar/enfermedades-cerebrales-medulares-y-nerviosas/trastornos-del-sistema-nervioso-aut%C3%B3nomo/introducci%C3%B3n-al-sistema-nervioso-aut%C3%B3nomo
[2] Torres, A., & Serrano, D. C. (2024, diciembre 10). Introducción al sistema nervioso periférico.
[3] (S/f). Ouraring.com. Recuperado el 2 de mayo de 2025, de https://support.ouraring.com/hc/es/articles/360025441974-Variabilidad-de-la-frecuencia-card%C3%ADaca
[4] (S/f-b). Udlap.mx. Recuperado el 2 de mayo de 2025, de http://catarina.udlap.mx/u_dl_a/tales/documentos/meie/osorio_s_a/capitulo2.pdf






