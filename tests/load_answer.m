function expected = load_answer(answer_file)
    answer_dir = get_answer_dir();            
    save_file = strcat(answer_dir, '/', answer_file);
    expected = load(save_file);
end 
