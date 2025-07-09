/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    // ¡Esta línea es CRUCIAL! Le dice a Tailwind dónde buscar tus clases.
    // Asegúrate de que esta ruta sea correcta para tus archivos JSX/TSX.
    "./src/**/*.{js,jsx,ts,tsx}", 
  ],
  theme: {
    extend: {},
  },
  plugins: [],
}
