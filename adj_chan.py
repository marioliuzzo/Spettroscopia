import total_counts

PATH = total_counts.PATH
NOME_SPETTRO = 'Ba133_1.txt' #modificare con il nome del file
PATH = total_counts.os.path.join(PATH, NOME_SPETTRO)

counts = total_counts.counts
channels = total_counts.channels
total_counts.plt.plot(channels, counts, marker = 'o')
NOME_SPETTRO = NOME_SPETTRO.replace('_2.txt', '')
NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')
total_counts.plt.title('Channels vs counts' + ' ' + NOME_SPETTRO)
total_counts.plt.xlabel('Channels')
total_counts.plt.ylabel('Counts')
total_counts.plt.minorticks_on()
total_counts.plt.show()
