import random 

player=int(input('Rock Paper Scissors 1.- ✊ Rock, 2.- ✋ Papper, 3.- ✌️ Scissors Selecciona un numero: '))
if player == 1:
	print('Tu escogiste: ✊ Rock')
elif player == 2:
	print('Tu escogiste: ✋ Papper')
elif player == 3:
	print('Tu escogiste: ✌️ Scissors')
else: 
	print('Seleccion equivocada')

	
CPU =random.randint(1,3)
if CPU==1:
	print("CPU escogió ✊ Rock")
elif CPU==2:
	print("CPU escogió ✋ Papper")
elif CPU==3:
	print("CPU escogió ✌️ Scissors")
else:
	print("Error de codigo")	

if player == CPU:
	print("Empate!")
elif player==1 and CPU==2:
	print("CPU gana!")
elif player==1 and CPU==3:
	print("player gana!")
elif player==2 and CPU==1:
	print("player gana!")
elif player==2 and CPU==3:
	print("CPU gana!")
elif player==3 and CPU==1:
	print("CPU gana!")
elif player==3 and CPU==2:
	print("player gana!")
else:
	print('Error de caso')

