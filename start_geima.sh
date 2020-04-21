#!/bin/bash
nohup /home/user/miniconda3/envs/django/bin/python /home/user/Web/imav2/manage.py runserver 172.16.0.55:8000 1>/home/user/Web/imav2/web.o 2>/home/user/Web/imav2/web.e &
