"""
Navier-Stokes Blowup Laboratory Launcher
Starts a local web server and opens the interactive 3D laboratory in your browser.

Usage:
  python3 scripts/run_lab.py
  python3 scripts/run_lab.py --port 8088
  python3 scripts/run_lab.py --test
"""

import os
import sys
import argparse
import webbrowser
import http.server
import socketserver
import threading
import time
import urllib.request


def find_free_port(starting_port=8080):
    import socket
    port = starting_port
    while port < starting_port + 100:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            res = sock.connect_ex(('127.0.0.1', port))
            if res != 0:
                return port
        port += 1
    return starting_port


class CustomHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        # Serve from project root directory
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        super().__init__(*args, directory=os.path.join(project_root, "web"), **kwargs)

    def end_headers(self):
        # Enable CORS and disable caching during interactive development
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate')
        super().end_headers()


def start_server(port: int, auto_open: bool = True):
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    web_dir = os.path.join(project_root, "web")

    if not os.path.exists(web_dir):
        print(f"Error: Web directory not found at {web_dir}")
        sys.exit(1)

    url = f"http://localhost:{port}"
    print("=" * 70)
    print("  Navier-Stokes Blowup Interactive 3D Laboratory")
    print(f"  Local URL: {url}")
    print(f"  Serving files from: {web_dir}")
    print("  Press Ctrl+C to stop the server.")
    print("=" * 70)

    # Allow socket reuse to prevent port-in-use errors upon restart
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("", port), CustomHandler) as httpd:
        if auto_open:
            threading.Timer(0.6, lambda: webbrowser.open(url)).start()
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nShutting down Navier-Stokes Laboratory server. Goodbye!")


def test_server(port: int):
    """Test mode: boots server, sends HTTP GET, verifies 200 OK, shuts down."""
    socketserver.TCPServer.allow_reuse_address = True
    httpd = socketserver.TCPServer(("", port), CustomHandler)
    server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    server_thread.start()
    time.sleep(0.5)

    test_url = f"http://localhost:{port}/index.html"
    try:
        req = urllib.request.urlopen(test_url, timeout=3)
        status = req.getcode()
        content = req.read().decode('utf-8')
        assert status == 200, f"Expected status 200, got {status}"
        assert "Navier–Stokes" in content, "Expected title in index.html content"
        print(f"[TEST PASSED] Server responded with status 200 OK from {test_url}")
    finally:
        httpd.shutdown()
        httpd.server_close()


def main():
    parser = argparse.ArgumentParser(description="Navier-Stokes Laboratory Web Server")
    parser.add_argument("--port", type=int, default=None, help="Port to listen on")
    parser.add_argument("--no-open", action="store_true", help="Do not open browser automatically")
    parser.add_argument("--file", action="store_true", help="Open index.html directly as a local file URL")
    parser.add_argument("--test", action="store_true", help="Run self-test and exit")
    args = parser.parse_args()

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    index_path = os.path.join(project_root, "web", "index.html")

    if args.file:
        print(f"Opening {index_path} directly in your browser...")
        webbrowser.open(f"file://{index_path}")
        return

    port = args.port if args.port is not None else find_free_port(8080)

    if args.test:
        test_server(port)
    else:
        start_server(port, auto_open=not args.no_open)


if __name__ == "__main__":
    main()
