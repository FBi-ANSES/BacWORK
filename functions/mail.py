import json
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

def get_projects(config):
	projects = []
	for sample in config["samples"]:
		project = sample["Project"]
		if project not in projects :
			projects.append(project)
	return ','.join(projects)

def get_samples(config):
	sample_list = []
	if config["Analysis"]["resume"] == "yes":
		json_files = resume_files_list(config)
		for json_file in json_files :
			if "summary.json" not in json_file:
				continue
			dico = json.load(open(json_file,'r'))
			sample_id = dico["SampleID"]
			quality = "PASS"
			for element in dico["Quality_filter"]:
				if dico["Quality_filter"][element] == "FALL":
					quality = "FALL"
					continue
			sample_list.append(f"- {sample_id}\t\t\t({quality})")		
	else :
		for sample in config["samples"]:
			sample_list.append(sample["SampleID"])
	return '\n'.join(sample_list)

def send_mail_onsuccess(config):
	if config["mail"]["send_mail"] == "yes" :
		
		server = config["mail"]["server"]
		port = config["mail"]["port"]
		mailfrom = config["mail"]["mailfrom"]
		rcptto = config["mail"]["recipients"]
		PROJECTS = get_projects(config)
		SAMPLES = get_samples(config)
		
		bilan_email = f"Bonjour,\n\nLes analyses pour le(s) projet(s) {PROJECTS} sont finies.\n"
		bilan_email = f"{bilan_email}\nLes echantillons analyses:\n{SAMPLES}"
		msg = MIMEMultipart("alternative")
		msg["Subject"] = f"[{PROJECTS}] analyse finies"
		msg["From"] = mailfrom
		msg["To"] = ",".join(rcptto)
		bilan_email = f"{bilan_email}\n\nBonne journee"
		text = bilan_email
		part1 = MIMEText(text, "plain")
		msg.attach(part1)

		try:
			smtp = smtplib.SMTP(server, port)     # Connexion au serveur
			smtp.sendmail(mailfrom, rcptto, msg.as_string())      # Envoi du mail
		except Exception as e: # Echec de l envoi
			print ('Failed to send mail.')
			print(str(e))
		else:
			print('Succeeded to send mail.')
		finally:
			if smtp != None:
				smtp.close()
				
def send_mail_onerror(config):
	if config["mail"]["send_mail"] == "yes" :
	
		server = config["mail"]["server"]
		port = config["mail"]["port"]
		mailfrom = config["mail"]["mailfrom"]
		rcptto = config["mail"]["recipients"]
		PROJECTS = get_projects(config)
		
		msg = MIMEMultipart("alternative")
		msg["Subject"] = f"[{PROJECTS}] soucis analyses"
		msg["From"] = mailfrom
		msg["To"] = ",".join(rcptto)
		bilan_email = f"Probleme analyse"
		text = bilan_email
		part1 = MIMEText(text, "plain")
		msg.attach(part1)
	 
		try:
			smtp = smtplib.SMTP(server, port)     # Connexion au serveur
			smtp.sendmail(mailfrom, rcptto, msg.as_string())      # Envoi du mail
		except Exception as e: # Echec de l envoi
			print ('Failed to send mail.')
			print(str(e))
		else:
			print('Succeeded to send mail.')
		finally:
			if smtp != None:
				smtp.close()
		print("Erreur snakemake")