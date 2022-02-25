function store_answer_file(answer_file, expected)

    answer_dir = get_answer_dir();

    save_file = strcat(answer_dir, answer_file);
    save(save_file, '-struct', 'expected')
end 
