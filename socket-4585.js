() => {
    const socket = new WebSocket('ws://emmago.hopto.org/programs/ws');
    socket.onopen = () => {
      console.log("el socket ha sido abierto")
    }
    socket.onerror = () => {
      console.log("hooooli")
    }
    return socket;
};