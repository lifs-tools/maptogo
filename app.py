from analysis import app
server = app.server

# Run the app
if __name__ == '__main__':
    app.run_server(host = "0.0.0.0", debug = True, port = 8040)
    #app.run_server(host = "0.0.0.0", debug = True, port = 80)
