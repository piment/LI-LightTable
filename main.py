#-*-coding:utf-8-*-
from modules.Adafruit_I2C import Adafruit_I2C

from modules.MCP23017 import MCP23017
from modules.lighttable import LightTable
import smbus
import time
import modules.opc

table = LightTable()
bus  = smbus.SMBus(0)

box = False

while True:
  # Attente de connexion du boitier avant demarrage du programme
  while not box:
      try:
          bus.read_byte(0x20)
          box = True
          print('Le boitier est connecte!')

      except:
          print("En attente de connexion du boitier de commande")
          print("Nouvel essaie dans 2s")
          time.sleep(2)

# Configuration du boitier de commande (MCP23017)
  try:
      table.client.put_pixels(table.AllOff)

      mcp = MCP23017(address=0x20, num_gpios=16, busnum=0)
      mcp.configSystemInterrupt(mcp.INTMIRRORON, mcp.INTPOLACTIVEHIGH)

      for i in range(0,8):
          mcp.pinMode(i, mcp.INPUT)
          mcp.pullUp(i, 1)
          mcp.configPinInterrupt(i, mcp.INTERRUPTON, mcp.INTERRUPTCOMPAREPREVIOUS)
      for i in range (8,16):
          mcp.pinMode(i, mcp.OUTPUT)
          mcp.output(i, 1)

      mcp.output(15, 0)
      mcp.output(8, 0)

      position = 8
      print('Le boitier est prÃªt')

  except:
      pass

  last_case = [0, 0, 0, 0]
  actual_cases = [0, 0, 0, 0]
  game = True
  while game:

    # On lit l'interrupt
      try:
          pin, value = mcp.readInterrupt()
      except:
          game = False
    # Si Pin a une valeur, alors un des boutons a ete presse
      if pin:

        # Si le Bouton Valider est presse©
          if pin == 7 and value == 0:

              for id, case in enumerate(actual_cases):
                  last_case[id] = actual_cases[id]

              print("Avant Lecture Case = ", actual_cases)
              my_inputs = ""
              for i in range(0, 5):
                  my_inputs += str(mcp.input(i))
              case = int(my_inputs, 2)
              actual_cases[position - 8] = case
              print("AprÃ©s lecture Case = ", actual_cases)
              print("Last case = ", last_case)
              for case in last_case:
                  if case:
                      for pos in table.get_case(case):
                          if pos:
                              table.lines[pos] = table.color_off
              table.client.put_pixels(table.AllOff)
              temp = 0   
              for index, case in enumerate(actual_cases):
                                     
                  if case <= 20:
                      for index2, pos in enumerate(table.get_case(actual_cases[index])):
                         print("case = " + str(case) + " and pos = " + str(pos))
                         if case > 0:
                             if pos == 0:
                                 temp = 1
                                 table.lines[pos] = table.color_green
                             else:
                                 table.lines[pos] = table.color_green
                         # if case == 0 and temp == 0:
                             #table.lines[pos] = table.color_off
                            
                         else:
                             if temp == 1:
                                 pass
                             else:
                                 table.lines[pos] = table.color_off
                             # temp = 1 #table.lines[0] = table.color_off
                  else:
                      table.client.put_pixels(table.Error)
                      time.sleep(0.5)
                      table.client.put_pixels(table.AllOff)

              table.client.put_pixels(table.lines)


        # Si le Bouton RAZ est presse
          if pin == 6 and value == 0:
              mcp.output(position, 1)
              for a in range(0, 2):
                  for b in range(8,12):
                      table.client.put_pixels(table.AllOn)
                      mcp.output(b, 0)
                      time.sleep(0.1)
                      table.client.put_pixels(table.AllOff)
                      mcp.output(b, 1)
                      time.sleep(0.1)
              position = 8
              mcp.output(position, 0)
              last_case = [0, 0, 0, 0]
              actual_cases = [0, 0, 0, 0]
              table.client.put_pixels(table.AllOff)
              table.table_raz()
              print("Last Case after RAZ = ",last_case)


            # TODO : Efface toute les cases presentes sur la table

          if pin == 5 and value == 0:
              mcp.output(position, 1)
              mcp.output(position + 1, 0)
              if position + 1 == 12:
                  mcp.output(position, 1)
                  position = 8
                  mcp.output(position, 0)
              else:
                  position += 1
            # TODO : Selectionne la position suivante dans le tableau

      mcp.clearInterrupts()
      time.sleep(1)
      """
    pin, value = mcp.readInterrupt()
    if pin != None:
        if pin == 0:
            if position:
                position = False
                mcp.output(8, 0)
            else:
                position = True
                mcp.output(8, 1)
            print('Changement de position')
        elif pin == 1:
            print('Remise Ã  zÃ©ro des position')
        elif pin == 2:
            print('Validation des coordonnÃ©es')

    time.sleep(2)
    """
      """  
    print("pin = {} / Value = {}".format(pin1, value1))
    if mcp.input(8) == 0:
        mcp.output(9,1)
    else:
        mcp.output(9, 0)
    time.sleep(1)	
    """
