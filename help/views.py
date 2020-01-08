from django.shortcuts import render

# Create your views here.

def index(request):
    return render(request, 'help.html', {'title': "系统生物学平台"})