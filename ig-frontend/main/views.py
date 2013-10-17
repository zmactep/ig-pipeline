from django.shortcuts import render


def about(request):
    return render(request, "main/about.html")


def contact(request):
    if not request.POST:
        return render(request, "main/contact.html")
    else:
        return render(request, "main/sent.html")
